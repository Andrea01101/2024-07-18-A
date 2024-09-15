import copy

import networkx as nx

from database.DAO import DAO
from model.gene import Gene


class Model:
    def __init__(self):
        self._graph = nx.DiGraph()
        self._localization_map = {}
        self._correlations_map = {}
        self.get_all_correlations()
        self._cammino_ottimo = []
        self._peso_ottimo = 0
        self._idMap = {}

    def get_all_correlations(self):
        DAO.get_all_correlations(self._correlations_map)

    def get_chromosomes(self):
        return DAO.get_all_chromosomes()

    def get_localizations(self):
        return DAO.get_all_localizations()

    def get_peso(self, g1: Gene, g2: Gene):
        return DAO.get_peso(g1, g2)

    def get_localization_gene(self, g: Gene):
        if g.GeneID in self._localization_map:
            return self._localization_map[g.GeneID]
        else:
            return DAO.get_localization_gene(g, self._localization_map)

    def build_graph(self, ch_min, ch_max):
        self._graph.clear()
        nodes = DAO.get_nodes(ch_min, ch_max)
        self._graph.add_nodes_from(nodes)
        for n in nodes:
            self._idMap[n.GeneID] = n

        for i in range(len(nodes) - 1):
            for j in range(i + 1, len(nodes)):
                if (self.get_localization_gene(nodes[i]) == self.get_localization_gene(nodes[j]) and
                        nodes[i].GeneID != nodes[j].GeneID and
                        (nodes[i].GeneID, nodes[j].GeneID) in self._correlations_map):
                    peso = self._correlations_map[(nodes[i].GeneID, nodes[j].GeneID)]
                    if nodes[i].Chromosome < nodes[j].Chromosome:
                        self._graph.add_edge(nodes[i], nodes[j], weight=peso)
                    elif nodes[i].Chromosome > nodes[j].Chromosome:
                        self._graph.add_edge(nodes[j], nodes[i], weight=peso)
                    else:
                        self._graph.add_edge(nodes[i], nodes[j], weight=peso)
                        self._graph.add_edge(nodes[j], nodes[i], weight=peso)




    def num_nodes(self):
        return len(self._graph.nodes)

    def nodes(self):
        return self._graph.nodes

    def num_edges(self):
        return len(self._graph.edges)

    def edges(self):
        return self._graph.edges

    def get_node_max_uscenti(self):
        sorted_nodes = sorted(self._graph.nodes(), key=lambda n: self._graph.out_degree(n), reverse=True)
        lista = []
        listaTupla = []
        for i in range(0, 5):
            lista.append(sorted_nodes[i])

        for n in lista:
            somma = 0
            for v in self._graph.neighbors(n):
                somma += self._graph[n][v]["weight"]
            listaTupla.append((n, somma))

        return listaTupla




        #result = []
        #for i in range(min(len(sorted_nodes), 5)):
         #   peso_tot = 0.0
          #  for e in self._graph.out_edges(sorted_nodes[i], data=True):
           #     peso_tot += float(e[2].get("weight"))
            #result.append((sorted_nodes[i], self._graph.out_degree(sorted_nodes[i]), peso_tot))
        #return result

    def get_connesse(self):
        return nx.weakly_connected_components(self._graph)

    def get_nodes_location(self, loc):
        res = []
        nodes = list(self._graph.nodes())
        for n in nodes:
            if self._localization_map[n.GeneID] == loc:
                res.append(n)
        return res

    def trova_cammino(self):
        self._cammino_ottimo = []
        self._peso_ottimo = 0
        for n in self._graph.nodes():
            nuovi_successori = self._calcola_successori_ammissibili(n, [n])
            self._ricorsione([n], nuovi_successori)
        return self._cammino_ottimo, self._peso_ottimo

    def _ricorsione(self, parziale, successori):
        # caso terminale
        if len(successori) == 0:
            if len(parziale) > len(self._cammino_ottimo):
                self._cammino_ottimo = copy.deepcopy(parziale)
                self._peso_ottimo = self._peso_cammino(self._cammino_ottimo)
            elif len(parziale) == len(self._cammino_ottimo) and self._peso_cammino(parziale) < self._peso_ottimo:
                self._cammino_ottimo = copy.deepcopy(parziale)
                self._peso_ottimo = self._peso_cammino(self._cammino_ottimo)
        # caso ricorsivo
        else:
            for n in successori:
                parziale.append(n)
                nuovi_successori = self._calcola_successori_ammissibili(n, parziale)
                self._ricorsione(parziale, nuovi_successori)
                parziale.pop()

    def _calcola_successori_ammissibili(self, n, parziale):
        last_essential = parziale[-1].Essential
        if len(parziale) == 1:
            nuovi_successori = [i for i in list(self._graph.successors(n)) if
                                i not in parziale and i.Essential != last_essential]
        else:
            last_peso = self._graph.get_edge_data(parziale[-2], parziale[-1])["weight"]
            nuovi_successori = [i for i in list(self._graph.successors(n)) if
                                i not in parziale and i.Essential != last_essential
                                and self._graph.get_edge_data(parziale[-1], i)["weight"] >= last_peso]
        return nuovi_successori

    def _peso_cammino(self, cammino):
        peso = 0
        if len(cammino) == 1:
            return peso
        for i in range(0, len(cammino) - 1):
            peso += self._graph.get_edge_data(cammino[i], cammino[i + 1])["weight"]
        return peso
