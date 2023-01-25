from oakvar import BasePostAggregator
from pathlib import Path
import sys
cur_path = str(Path(__file__).parent)
sys.path.append(cur_path)
import sqlite3
import longevitymap_ref_homo
import json

MULTIPLE_CONST = "multiple"
CONFLICTED_CONST = "conflicted"
CONFLICTED_INDEX = -1

class CravatPostAggregator (BasePostAggregator):
    sql_insert:str = """ INSERT INTO longevitymap (
                weight,
                weightcolor,
                population,
                snp,
                gene,
                conflicted_rows,
                description,
                coding,
                ref,
                alt,
                cdnachange,
                deseases,
                zegot,
                alelfreq,
                nucleotides,
                priority,
                ncbidesc         
            ) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?) """

    ref_homo:longevitymap_ref_homo.RefHomoEdgecases = longevitymap_ref_homo.RefHomoEdgecases()

    def check(self):
        return True


    def get_nucleotides(self, ref:str, alt:str, zygocity:str) -> (str, set[str]):
        if zygocity == 'hom':
            return alt+"/"+alt, {alt, alt}
        return alt+"/"+ref, {alt, ref}


    def get_color(self, w: float, scale: float = 1.5) -> str:
        w = float(w)
        if w < 0:
            w = w * -1
            w = 1 - w * scale
            if w < 0:
                w = 0
            color: str = format(int(w * 255), 'x')
            if len(color) == 1:
                color = "0" + color
            color = "ff" + color + color
        else:
            w = 1 - w * scale
            if w < 0:
                w = 0
            color = format(int(w * 255), 'x')
            if len(color) == 1:
                color = "0" + color
            color = color + "ff" + color

        return color


    def cleanup (self):
        if self.result_cursor is not None:
            self.result_cursor.close()
        if self.result_conn is not None:
            self.result_conn.commit()
            self.result_conn.close()
        return


    def setup (self):
        self.result_path:Path = Path(self.output_dir, self.run_name + "_longevity.sqlite")
        self.result_conn:sqlite3.Connection = sqlite3.connect(self.result_path)
        self.result_cursor:sqlite3.Cursor = self.result_conn.cursor()


        sql_create:str = """ CREATE TABLE IF NOT EXISTS longevitymap (
            id integer NOT NULL PRIMARY KEY,
            weight float,
            weightcolor float,
            population text,
            snp text,
            gene text,
            conflicted_rows text,
            description text,
            coding text,
            ref text,
            alt text,
            cdnachange text,
            deseases text,
            zegot text,
            alelfreq text,
            nucleotides text,
            priority float,
            ncbidesc text          
            )"""
        self.result_cursor.execute(sql_create)
        self.result_cursor.execute("DELETE FROM longevitymap;")
        self.result_conn.commit()

        cur_path:str = str(Path(__file__).parent)
        self.data_conn:sqlite3.Connection = sqlite3.connect(Path(cur_path, "data", "longevitymap.sqlite"))
        self.data_cursor:sqlite3.Cursor = self.data_conn.cursor()
        # self.ref_homo.init(self, self.longevity_cursor, self.sql_insert)
        self.ref_homo.setup(self, self.result_cursor, self.data_cursor, self.sql_insert)


    def merge_records(self, row:tuple, record:list) -> list:
        need_info:bool = False
        if record is None:
            record = list(row)
            record[6] = []
            record[7] = []

        record[7].append({"pubmedid":row[5], "study_design":row[6], "conclusions":row[7]})

        if len(record) < 8:
            print("Error les than 7 len ----------------------------------------")
        if len(record[7]) == 0:
            print("Error 7 is empty----------------------------------------------")

        #id
        if record[0] != row[0]:
            record[0] = MULTIPLE_CONST
            need_info = True

        #association
        if record[1] != row[1]:
            record[1] = CONFLICTED_CONST
            need_info = True

        #population
        if record[2] != row[2]:
            record[2] = MULTIPLE_CONST
            need_info = True
        #identifier
        if record[3] != row[3]:
            record[3] = CONFLICTED_INDEX
            need_info = True

        #symbol (GENE)
        if record[4] != row[4]:
            record[4] = MULTIPLE_CONST
            need_info = True

        #quickpubmed
        if record[5] != row[5]:
            record[5] = CONFLICTED_INDEX
            need_info = True

        if need_info:
            # print("Info--------------------------------------")
            record[6].append({"id":row[0], "association":row[1], "population":row[2], "identifier":row[3], "gene":row[4], "pubmedid":row[5]})

        return record


    # 'longevitydb_id': str(record[0]),
    # 'association': str(record[1]),
    # 'population': str(record[2]),
    # 'rsid': str(record[3]),
    # 'genes': str(record[4]),
    # 'pmid': str(record[5]),
    # 'info': str(record[6]),
    # 'description': str(record[7]),
    # 'allele': str(allel_row[0]),
    # 'state': str(allel_row[1]),
    # 'zygosity': str(allel_row[2]),
    # 'weight': str(allel_row[3]),
    # 'priority': str(priority)

    def annotate (self, input_data):
        rsid:str = str(input_data['dbsnp__rsid'])
        if rsid == '':
            return

        if not rsid.startswith('rs'):
            rsid = "rs" + rsid
        query:str = 'SELECT variant.id, association, population.name, identifier, symbol, quickpubmed, study_design, conclusions ' \
                'FROM variant, population, gene, allele_weights WHERE  ' \
                'variant.identifier = "{rsid}" AND variant.population_id = population.id AND variant.gene_id = gene.id AND ' \
                'allele_weights.rsid = variant.identifier AND allele_weights.allele = "{alt}" GROUP BY variant.id'.format(
            rsid=rsid, alt=input_data['base__alt_base'])

        self.data_cursor.execute(query)
        rows:tuple = self.data_cursor.fetchall()

        if len(rows) == 0:
            return None

        record:list = None
        for row in rows:
            record = self.merge_records(row, record)

        zygot:str = input_data['vcfinfo__zygosity']
        if zygot is None or zygot == "":
            zygot = "het"

        alt:str = input_data['base__alt_base']
        ref:str = input_data['base__ref_base']

        query2:str = f"SELECT weight, priority FROM allele_weights WHERE rsid = '{rsid}' AND zygosity = '{zygot}' AND allele = '{alt}'"
        self.data_cursor.execute(query2)
        rows2:tuple = self.data_cursor.fetchall()
        if len(rows2) == 0:
            return
        allel_row:tuple = rows2[0]
        w = allel_row[0]
        priority:str = allel_row[1]

        if len(rows2) > 1:
            print("Worning unexpected number of rows in allel_row in longevitymap postagregator!!!____________________________________________")

        if record[1] != "significant":
            return

        self.ref_homo.process_row(input_data)
        nuq, nuq_set = self.get_nucleotides(ref, alt, zygot)

        if w == 0:
            return

        color:str = self.get_color(w, 1.5)
        # temp = self._createSubTable(record[6])
        # temp += record[7].replace("____", "<br/>").replace("__", " ")

        task:tuple = (w, color, record[2], rsid, record[4], json.dumps(record[6]), json.dumps(record[7]),
                input_data['base__coding'], ref, alt, input_data['base__cchange'], input_data['clinvar__disease_names'],
                zygot, input_data['gnomad__af'], nuq, priority, input_data['ncbigene__ncbi_desc'])

        self.result_cursor.execute(self.sql_insert, task)
        return {"col1":""}


    def postprocess(self):
        self.ref_homo.end()
