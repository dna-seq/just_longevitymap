from oakvar import BasePostAggregator
from pathlib import Path
import sys
cur_path = str(Path(__file__).parent)
sys.path.append(cur_path)
import sqlite3
import longevitymap_ref_homo
import json
import requests
import csv

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
                ncbidesc,
                category_name
            ) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?) """

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
        self.categories = set(self.categories)
        self.get_llm_answer()
        if self.result_cursor is not None:
            self.result_cursor.close()
        if self.result_conn is not None:
            self.result_conn.commit()
            self.result_conn.close()
        return


    def setup (self):
        self.categories = []
        self.result_path:Path = Path(self.output_dir, self.run_name + "_longevity.sqlite")
        self.result_conn:sqlite3.Connection = sqlite3.connect(self.result_path)
        self.result_cursor:sqlite3.Cursor = self.result_conn.cursor()

        self.json_path = Path(self.output_dir, "output.json")
        if self.json_path.is_file():
            pass
        else:
            with open(self.json_path, "w") as json:
                json.write("[{}]")


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
            ncbidesc text,
            category_name test
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

        if len(record) < 9:
            print("Error les than 8 len ----------------------------------------")
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

        #category_name
        if record[8] != row[8]:
            record[8] = CONFLICTED_INDEX
            need_info = True

        if need_info:
            # print("Info--------------------------------------")
            record[6].append({"id":row[0], "association":row[1], "population":row[2], "identifier":row[3], "gene":row[4], "pubmedid":row[5], "category_name":row[8]})

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


    def get_llm_answer(self):
        for category in self.categories:
            tsv_path = Path(self.output_dir, category + ".tsv")
            cursor_to_sqlite = self.result_cursor
            cursor_to_sqlite.execute(f"SELECT * FROM longevitymap WHERE category_name = '{category}';")
            column_names = list(map(lambda x: x[0], cursor_to_sqlite.description))
            rows = cursor_to_sqlite.fetchall()

            with open(tsv_path, 'w', encoding="utf-8") as tsv_file:
                writer = csv.writer(tsv_file, delimiter='\t')
                writer.writerow(column_names)
                for row in rows:
                    writer.writerow(list(row))

            json_path = self.json_path
            tsv_path = Path(self.output_dir, category + ".tsv")

            url = "http://agingkills.eu:8088/v1/chat/completions"

            prompt = """
                        You will be provided a piece of personalized genomic report of some person.
                        You have to analyze information and based on it provide recommendations that can be made based on the information.
                        This recommendations should contain the following information:
                        Identification of Crucial Genotypes: Highlighting specific genotypes that are significant for the individual's health,
                        such as those related to common genetic disorders, responses to medications, or predispositions to certain conditions.
                        Devide those Crucial Genotypes into two groups: 'Beneficial Genotypes' and 'Detrimental Genotypes'.
                        Interactions Between Genotypes: Analyzing how different genotypes might interact with each other, either amplifying or
                        mitigating risks. For example, some genetic variants might work together to increase the risk of a particular disease.
                        Risk Assessment: Providing information on the risks associated with certain genotypes, such as the likelihood of
                        developing specific diseases or conditions, always devide those risks into two groups: longevity related risks and general risks.
                        Risk Reduction Strategies: Offering recommendations for lifestyle changes,
                        dietary adjustments, or medical interventions that could help mitigate the risks associated with certain genetic profiles.
                        Additional Insights: Any other relevant information that could be
                        useful for the individual based on their genomic profile.
                        Always make information clear, talk only about genes and genotypes that were provided to you, It should be clear, that you provide
                        the summarization of the information of the genotypes and genes that were provided to you in the report.
                        In your answer do not write the header of the whole text, for example, do not write "Personalized Genomic Report Analysis".
                        Do not ask follow-up questions. Do not ask if the user understood the provided information or not.
                        Also you should write your answer in HTML formatting. For small headers like "Risk Assessment", "Identification of Crucial Genotypes",
                        "Interactions Between Genotypes", "Risk Reduction Strategies" and "Additional Insights" use tag <h3> with class="small-headers".
                        When you want to provide a link you should use html tag <a> with "href" attribute where you place the link.
                        If you want to make some of the text emphasized use <b> tag. Do not forget to divid information in different paragraphs so it looks
                        better in HTML page. Never forget to close tags. Your answer will be inserted in already existing HTML page. So make sure to format
                        it accordingly for this purpose. Everything should be wrapped in tags. If you don't know what tag to use for the main text, use <p> tag.
                        Please use the following "summ" class for the div tag that you use to wrap the whole your answer.
                        Here is the report information in tsv format style: '''
                    """
            with open(tsv_path, "r") as tsv_file:
                tsv_info = tsv_file.read()

            text_for_request = prompt + str(tsv_info)

            json_api_openai = {'model': 'gpt-4o',
                            'messages': [{'role': 'system', 'content': ''},
                                            {'role': 'user', 'content': [{'type': 'text', 'text': text_for_request}]}
                                            ],
                                'stream': False,
                                'max_tokens': 10000,
                                'stop': ['[DONE]'],
                                'temperature': 0}

            json_for_mistral = {'model': 'mistral-large-latest',
                                'messages': [{'role': 'system', 'content': ''},
                                            {'role': 'user', 'content': [{'type': 'text', 'text': text_for_request}]}],
                                            'stream': True, 'max_tokens': 10000,
                                            'stop': ['[DONE]'],
                                            'temperature': 0}
            try:
                answer = requests.post(url, json=json_api_openai)
            except:
                print(f"Can't make request to {url}")

            answer = answer.json()
            answer = answer["choices"][0]["message"]["content"]

            answer = {category: str(answer)}

            with open (json_path, mode="r+") as json_file:
                json_file.seek(0,2)
                position = json_file.tell() -1
                json_file.seek(position)
                json_file.write( ",{}]".format(json.dumps(answer)))
        return


    def annotate (self, input_data):
        rsid:str = str(input_data['dbsnp__rsid'])
        if rsid == '':
            return

        if not rsid.startswith('rs'):
            rsid = "rs" + rsid
        query:str = 'SELECT variant.id, association, population.name, identifier, symbol, quickpubmed, study_design, conclusions, categories.name ' \
                'FROM variant, population, gene, allele_weights, categories WHERE  ' \
                'variant.identifier = "{rsid}" AND variant.population_id = population.id AND variant.gene_id = gene.id AND ' \
                'allele_weights.rsid = variant.identifier AND allele_weights.allele = "{alt}" AND allele_weights.category_id=categories.id GROUP BY variant.id'.format(
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
            print("Warning unexpected number of rows in allel_row in longevitymap postagregator!!!____________________________________________")

        if record[1] != "significant":
            return
        
        self.ref_homo.process_row(input_data)
        nuq, nuq_set = self.get_nucleotides(ref, alt, zygot)

        if w == 0:
            return
        if len(str(w))>5:
            w = round(w, 3)
        color:str = self.get_color(w, 1.5)
        # temp = self._createSubTable(record[6])
        # temp += record[7].replace("____", "<br/>").replace("__", " ")

        self.categories.append(record[8])
        task:tuple = (w, color, record[2], rsid, record[4], json.dumps(record[6]), json.dumps(record[7]),
                input_data['base__coding'], ref, alt, input_data['base__cchange'], input_data['clinvar__disease_names'],
                zygot, input_data['gnomad__af'], nuq, priority, input_data['ncbigene__ncbi_desc'], record[8])

        self.result_cursor.execute(self.sql_insert, task)
        return {"col1":""}


    def postprocess(self):
        self.ref_homo.end()


