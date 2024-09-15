import flet as ft
from UI.view import View
from model.modello import Model


class Controller:
    def __init__(self, view: View, model: Model):
        # the view, with the graphical elements of the UI
        self._view = view
        # the model, which implements the logic of the program and holds the data
        self._model = model

    def handle_graph(self, e):
        #read user input
        if self._view.dd_min_ch.value is None:
            self._view.create_alert("Selezionare un valore per Chromosoma min")
            return
        else:
            ch_min = int(self._view.dd_min_ch.value)
        if self._view.dd_max_ch.value is None:
            self._view.create_alert("Selezionare un valore per Chromosoma max")
            return
        else:
            ch_max = int(self._view.dd_max_ch.value)
        if ch_min > ch_max:
            self._view.create_alert("Attenzione: deve essere Chromosoma min <= Chromosoma max")
            return

        # crea grafo e stampa info del grafo
        self._model.build_graph(ch_min, ch_max)
        self._view.txt_result1.controls.clear()
        self._view.txt_result1.controls.append(ft.Text(f"Creato grafo con {self._model.num_nodes()} nodi"
                                                       f" e {self._model.num_edges()} archi"))

        sorted_nodes = self._model.get_node_max_uscenti()
        for v, pesoTot in sorted_nodes:
            self._view.txt_result1.controls.append(ft.Text(f"{v}| weight: {pesoTot}"))

        self._view.btn_dettagli.disabled = False
        self._view.btn_path.disabled = False
        self._view.update_page()

    def handle_dettagli(self, e):
        loc = self._view.dd_localization.value
        nodes = self._model.get_nodes_location(loc)
        self._view.txt_result1.controls.clear()
        self._view.txt_result1.controls.append(ft.Text(f"\n Ci sono {len(nodes)} geni con location {loc}:"))
        for n in nodes:
            self._view.txt_result1.controls.append(ft.Text(f"{n}"))
        self._view.update_page()


    def handle_path(self, e):
        self._view.txt_result2.controls.clear()
        cammino_ottimo, costo_ottimo = self._model.trova_cammino()
        self._view.txt_result2.controls.append(ft.Text(f"Trovato un cammino ottimo di lunghezza {len(cammino_ottimo)}"))
        self._view.txt_result2.controls.append(ft.Text(f"Il costo del cammino ottimo è {costo_ottimo}"))
        self._view.txt_result2.controls.append(ft.Text(f"Le fermate del cammino sono:"))
        for n in cammino_ottimo:
            self._view.txt_result2.controls.append(ft.Text(f"{n}"))
        self._view.update_page()

    def fill_dd_ch(self):
        values = self._model.get_chromosomes()
        for ch in values:
            self._view.dd_min_ch.options.append(ft.dropdown.Option(key=ch, text=ch))
            self._view.dd_max_ch.options.append(ft.dropdown.Option(key=ch, text=ch))

        #self._view.dd_min_ch.options = list(map(lambda x: ft.dropdown.Option(key=x, text=x), values))
        #self._view.dd_max_ch.options = list(map(lambda x: ft.dropdown.Option(key=x, text=x), values))
        self._view.update_page()

    def fill_dd_localization(self):
        values = self._model.get_localizations()
        self._view.dd_localization.options = list(map(lambda x: ft.dropdown.Option(key=x, text=x), values))
