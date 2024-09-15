from database.DB_connect import DBConnect
from model.gene import Gene
from model.interaction import Interaction


class DAO():

    @staticmethod
    def get_all_genes():
        cnx = DBConnect.get_connection()
        result = []
        if cnx is None:
            print("Connessione fallita")
        else:
            cursor = cnx.cursor(dictionary=True)
            query = """SELECT * 
                    FROM genes"""
            cursor.execute(query)

            for row in cursor:
                result.append(Gene(**row))

            cursor.close()
            cnx.close()
        return result

    @staticmethod
    def get_all_interactions():
        cnx = DBConnect.get_connection()
        result = []
        if cnx is None:
            print("Connessione fallita")
        else:
            cursor = cnx.cursor(dictionary=True)
            query = """SELECT * 
                       FROM interactions"""
            cursor.execute(query)

            for row in cursor:
                result.append(Interaction(**row))

            cursor.close()
            cnx.close()
        return result

    @staticmethod
    def get_all_chromosomes():
        cnx = DBConnect.get_connection()
        result = []
        if cnx is None:
            print("Connessione fallita")
        else:
            cursor = cnx.cursor(dictionary=True)
            query = """SELECT DISTINCT g.Chromosome  
                        FROM genes g
                        ORDER BY g.Chromosome ASC"""
            cursor.execute(query)

            for row in cursor:
                result.append(row["Chromosome"])

            cursor.close()
            cnx.close()
        return result

    @staticmethod
    def get_all_localizations():
        cnx = DBConnect.get_connection()
        result = []
        if cnx is None:
            print("Connessione fallita")
        else:
            cursor = cnx.cursor(dictionary=True)
            query = """SELECT DISTINCT c.Localization  
                        FROM classification c 
                        ORDER BY c.Localization  ASC"""
            cursor.execute(query)

            for row in cursor:
                result.append(row["Localization"])

            cursor.close()
            cnx.close()
        return result

    @staticmethod
    def get_nodes(ch_min:int, ch_max:int):
        cnx = DBConnect.get_connection()
        result = []
        if cnx is None:
            print("Connessione fallita")
        else:
            cursor = cnx.cursor(dictionary=True)
            query = """SELECT *  
                        FROM genes g
                        WHERE g.Chromosome >= %s AND g.Chromosome <= %s"""
            cursor.execute(query, (ch_min, ch_max))

            for row in cursor:
                result.append(Gene(**row))

            cursor.close()
            cnx.close()
        return result

    @staticmethod
    def get_localization_gene(g: Gene, localization_map: dict):
        cnx = DBConnect.get_connection()
        result = None
        if cnx is None:
            print("Connessione fallita")
        else:
            cursor = cnx.cursor(dictionary=True)
            query = """SELECT c.Localization  
                        FROM classification c 
                        WHERE c.GeneID = %s """
            cursor.execute(query, (g.GeneID,))

            rows = cursor.fetchall()
            localization_map[g.GeneID] = rows[0]["Localization"]
            result = rows[0]["Localization"]

            cursor.close()
            cnx.close()
        return result

    @staticmethod
    def get_all_correlations(correlations_map: dict):
        cnx = DBConnect.get_connection()
        if cnx is None:
            print("Connessione fallita")
        else:
            cursor = cnx.cursor(dictionary=True)
            query = """SELECT *  
                        FROM interactions i """
            cursor.execute(query)

            for row in cursor:
                id1 = row["GeneID1"]
                id2 = row["GeneID2"]
                corr = row["Expression_Corr"]
                correlations_map[(id1, id2)] = corr
                correlations_map[(id2, id1)] = corr
            cursor.close()
            cnx.close()

    @staticmethod
    def get_peso(g1: Gene, g2: Gene):
        cnx = DBConnect.get_connection()
        result = None
        if cnx is None:
            print("Connessione fallita")
        else:
            cursor = cnx.cursor(dictionary=True)
            query = """SELECT i.Expression_Corr as corr 
                        FROM interactions i 
                        WHERE i.GeneID1 = %s and i.GeneID2 = %s """
            cursor.execute(query, (g1.GeneID, g2.GeneID))

            rows = cursor.fetchall()

            if len(rows) > 0:
                result = rows[0]

            cursor.close()
            cnx.close()
        return result

    @staticmethod
    def get_archi(ch_min: int, ch_max: int):
        cnx = DBConnect.get_connection()
        result = []
        if cnx is None:
            print("Connessione fallita")
        else:
            cursor = cnx.cursor(dictionary=True)
            query = """SELECT DISTINCTROW g1.GeneID as g1id, g1.`Function` , g2.GeneID as g2id, g2.`Function` , g1.Chromosome as ch1, g2.Chromosome as ch2, c1.Localization, c2.Localization, i.Expression_Corr as peso
                        FROM genes g1, genes g2, classification c1, classification c2, interactions i
                        WHERE g1.GeneID = c1.GeneID
                        and g2.GeneID = c2.GeneID
                        and c1.Localization = c2.Localization
                        and g1.GeneID <> g2.GeneID
                        and g1.Chromosome >= %s
                        and g1.Chromosome <= %s
                        and g2.Chromosome >= %s
                        and g2.Chromosome <= %s
                        and ((i.GeneID1=g1.GeneID AND i.GeneID2=g2.GeneID) or (i.GeneID1=g2.GeneID AND i.GeneID2=g1.GeneID))
                        having g1.Chromosome <= g2.Chromosome """
            cursor.execute(query, (ch_min, ch_max, ch_min, ch_max))

            for row in cursor:
                result.append((row["g1id"], row["g2id"], row["peso"]))

            cursor.close()
            cnx.close()
        return result
    @staticmethod
    def getAllEdges(idMap):
        conn = DBConnect.get_connection()

        result = []

        cursor = conn.cursor(dictionary=True)
        query = """select i.GeneID1 as g1, i.GeneID2 as g2, i.Expression_Corr as peso
                   from interactions i, classification c, classification c2
                   where i.GeneID1 != i.GeneID2 and i.GeneID1 = c.GeneId
                   and i.GeneID2 = c2.GeneID and c.Localization = c2.Localization"""

        cursor.execute(query)

        for row in cursor:
            if row["g1"] in idMap and row["g2"] in idMap:  # solo se esistono appendo
                result.append((idMap[row["g1"]], idMap[row["g2"]], row["peso"]))

        cursor.close()
        conn.close()
        return result

    @staticmethod
    def get_nodes2(ch_min: int, ch_max: int):
        cnx = DBConnect.get_connection()
        result = []
        if cnx is None:
            print("Connessione fallita")
        else:
            cursor = cnx.cursor(dictionary=True)
            query = """SELECT *  
                        FROM genes g
                        WHERE g.Chromosome >= %s AND g.Chromosome <= %s"""
            cursor.execute(query, (ch_min, ch_max))

            for row in cursor:
                result.append(Gene(**row))

            cursor.close()
            cnx.close()
        return result  # Non restituire una tupla, solo la lista di nodi

    @staticmethod
    def getAllEdges2(idMap):
        conn = DBConnect.get_connection()

        result = []

        cursor = conn.cursor(dictionary=True)
        query = """select i.GeneID1 as g1, i.GeneID2 as g2, i.Expression_Corr as peso
                   from interactions i, classification c, classification c2
                   where i.GeneID1 != i.GeneID2 and i.GeneID1 = c.GeneId
                   and i.GeneID2 = c2.GeneID and c.Localization = c2.Localization"""

        cursor.execute(query)

        for row in cursor:
            g1 = row["g1"]
            g2 = row["g2"]

            if g1 not in idMap:
                print(f"g1 {g1} non trovato in idMap")
            if g2 not in idMap:
                print(f"g2 {g2} non trovato in idMap")

            if g1 in idMap and g2 in idMap:  # solo se esistono appendo
                result.append((idMap[g1], idMap[g2], row["peso"]))

        cursor.close()
        conn.close()
        return result  # Anche qui restituiamo solo la lista di archi


if __name__ == "__main__":
    pass
    # print(DAO.get_localization_gene("G234347", {}))
