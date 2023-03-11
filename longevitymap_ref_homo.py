import sqlite3
from sqlite3 import Error
from pathlib import Path
import json

ALLELE = "allele"
EXIST = "exist"
WEIGHT = "weight"

MULTIPLE_CONST = "multiple"
CONFLICTED_CONST = "conflicted"
CONFLICTED_INDEX = -1

class RefHomoEdgecases:
    _is_active:bool = True
    sql_ref_homozygot:str = """SELECT rsid, allele, weight FROM allele_weights WHERE state = 'ref' AND zygosity = 'hom'"""
    ref_homo_map:dict[str, dict] = {}


    def setup(self, parent, result_cursor:sqlite3.Cursor, data_cursor:sqlite3.Cursor, sql_insert:str) -> None:
        self.result_cursor:sqlite3.Cursor = result_cursor
        self.data_cursor:sqlite3.Cursor = data_cursor
        self.sql_insert: str = sql_insert
        self.parent = parent

        try:
            self.data_cursor.execute(self.sql_ref_homozygot)
            ref_homozygots:tuple = self.data_cursor.fetchall()
            for row in ref_homozygots:
                self.ref_homo_map[row[0]] = {ALLELE:row[1], WEIGHT:row[2], EXIST:True}

        except Error as e:
            print(e)


    def setActive(self, active:bool) -> None:
        self._is_active:bool = active


    def process_record(self, rsid, allele, w) -> None:
        if not self._is_active:
            return
        query:str = 'SELECT variant.id, association, population.name, identifier, symbol, quickpubmed, study_design, conclusions, category_name ' \
                'FROM variant, population, gene, allele_weights, snps, categories WHERE  ' \
                'variant.identifier = "{rsid}" AND snps.rsid="{rsid}" AND variant.population_id = population.id AND variant.gene_id = gene.id AND ' \
                'allele_weights.rsid = variant.identifier AND allele_weights.allele = "{alt}" AND' \
                ' allele_weights.state = "ref" AND allele_weights.zygosity = "hom" AND category=category_id' \
                ' GROUP BY variant.id'.format(
            rsid=rsid, alt=allele)

        self.data_cursor.execute(query)
        rows:tuple = self.data_cursor.fetchall()

        if len(rows) == 0:
            return

        record:list = None
        for row in rows:
            record = self.parent.merge_records(row, record)

        alt:str = allele
        ref:str = allele
        zygot:str = "hom"
        nuq:str = allele + "/" + allele
        color:str = self.parent.get_color(w, 1.5)

        task:tuple = (w, color, record[2], rsid, record[4], json.dumps(record[6]), json.dumps(record[7]), "", ref, alt, "", "", zygot, "", nuq, "0", "", record[8])

        self.longevity_cursor.execute(self.sql_insert, task)



    def process_row(self, row:list):
        if not self._is_active:
            return
        if len(self.ref_homo_map) == 0:
            return
        rsid:str = str(row['dbsnp__rsid'])
        if rsid == '':
            return
        if not rsid.startswith('rs'):
            rsid = "rs"+rsid
        item:dict = self.ref_homo_map.get(rsid)
        if item:
            self.ref_homo_map[rsid][EXIST] = False


    def end(self):
        if not self._is_active:
            return
        for rsid in self.ref_homo_map:
            if self.ref_homo_map[rsid][EXIST]:
                self.process_record(rsid, self.ref_homo_map[rsid][ALLELE], self.ref_homo_map[rsid][WEIGHT])