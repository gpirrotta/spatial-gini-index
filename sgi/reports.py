import numpy as np
import pandas as pd


class Report:

    def __init__(self, sgi_engine):
        self._sgi_engine = sgi_engine
        self._report = None
        self._report_pvalue = None
        self._metric = sgi_engine.metric
        
    

    def _build_report(self):
        if self._sgi_engine.metric == 'SGI':
            self._build_report_sgi()
        elif self._sgi_engine.metric =='SGD':
            self._build_report_sgd()

    def _build_report_sgi(self):
        
        schema = {
            "Spatial Lag": int,
            "h-distance": float,
            "Links": int,
            "Contiguity Variability": float
        }

        number_lags = np.arange(1,len(self._sgi_engine.h_distances)+1)

        data = pd.DataFrame(
            np.c_[
                number_lags,
                self._sgi_engine.h_distances,
                self._sgi_engine.links,
                self._sgi_engine.variabilities
            ],
            columns=schema.keys(),
        ).astype(schema)

        self._set_links_name(data)

        data["% Links"] = data["Links"] / data["Links"].sum(
                axis=0
            )
        data["% Contiguity Variability"] = data[
            "Contiguity Variability"
        ] / data["Contiguity Variability"].sum(axis=0)
        data["Cumulative Links"] = data["% Links"].cumsum()
        data["Cumulative Contiguity Variability"] = data[
            "% Contiguity Variability"
        ].cumsum()
        data["tg (rad)"] = np.arctan(
            data["Cumulative Contiguity Variability"]
            / data["Cumulative Links"]
        )
        data["tg (degree)"] = (data["tg (rad)"] * 360) / (2 * np.pi)

        data["SGI"] = ""
        data["Time in seconds"] = ""
        
        data.iat[
            0, data.columns.get_loc("SGI")
        ] = self._sgi_engine._sgi
        data.iat[
            0, data.columns.get_loc("Time in seconds")
        ] = self._sgi_engine._time_sgi

        data.fillna("", inplace=True)
        self._report= data


    def _build_report_sgd(self):    

        schema = {
            "Spatial Lag": int,
            "h-distance": float,
            "Links": int,
            "Contiguity Variability": float,
        }

        number_lags = np.arange(1,len(self._sgi_engine._h_distances)+1)

        data = pd.DataFrame(
            np.c_[
                number_lags,
                self._sgi_engine.h_distances,
                self._sgi_engine.links,
                self._sgi_engine.variabilities,
            ],
            columns=schema.keys(),
        ).astype(schema)

        self._set_links_name(data)

        data["% Links"] = data["Links"] / data["Links"].sum(
            axis=0
        )
        data["% Contiguity Variability"] = data[
            "Contiguity Variability"
        ] / data["Contiguity Variability"].sum(axis=0)

        data["Cumulative Links"] = data["% Links"].cumsum()
        data["Cumulative Contiguity Variability"] = data[
            "% Contiguity Variability"
        ].cumsum()        

        #Non contiguity
        data["Non Contiguity Links"] = (
                data["Links"].sum(axis=0) - data["Links"]
            )
        data["Non Contiguity Variability"] = (
                data["Contiguity Variability"].sum()
                - data["Contiguity Variability"]
            )

        data["% Non Contiguity Links"] = data[
                "Non Contiguity Links"
            ] / data["Non Contiguity Links"].sum(axis=0)
        
        data["% Non Contiguity Variability"] = data[
            "Non Contiguity Variability"
        ] / data["Non Contiguity Variability"].sum(axis=0)
        data["Cumulative Non Contiguity Links"] = data[
            "% Non Contiguity Links"
        ].cumsum()
        data["Cumulative Non Contiguity Variability"] = data[
            "% Non Contiguity Variability"
        ].cumsum()

        data["G (C)"] = data["Contiguity Variability"] / (
            2 * np.mean(self._sgi_engine._y) * self._sgi_engine._y.shape[0] ** 2
        )
        data["G (NC)"] = data["Non Contiguity Variability"] / (
            2 * np.mean(self._sgi_engine._y) * self._sgi_engine._y.shape[0] ** 2
        )
        data["G (ASP)"] = data["Contiguity Variability"].sum(axis=0) / (
            2 * np.mean(self._sgi_engine._y) * self._sgi_engine._y.shape[0] ** 2
        )

        data.fillna("", inplace=True)

        self._report = data


    def _set_links_name(self, data):
        if self._sgi_engine._step == "max-min" and self._sgi_engine._labels is not None:

            max_indexes = [np.where(self._sgi_engine._distance_matrix == max_threshold)[0] for max_threshold in self._sgi_engine._h_distances ]

            max_indexes_sr = pd.Series(max_indexes)

            tmp = pd.DataFrame()
            tmp['Place A'] = max_indexes_sr.map(lambda x: self._sgi_engine._labels[x[0]])
            tmp['Place B'] = max_indexes_sr.map(lambda x: self._sgi_engine._labels[x[1]])
            
            data['h-distance (links_name)'] = tmp['Place A'] +'-'+ tmp['Place B']



    def _build_report_pvalue(self):
        
        schema = {
            "Permutation": int,
            "SGI": float
        }


        data = pd.DataFrame(
            np.c_[
                np.arange(1,self._sgi_engine._permutations+1),
                self._sgi_engine._sgis            ],
            columns=schema.keys(),
        ).astype(schema)

        data["pvalue_sgi_z"] = ""
        data["zcal"] = ""

        data["pvalue_sgi_mc"] = ""


        data["Time (seconds) p-value SGI"] = ""
        

        data.iat[
            0, data.columns.get_loc("pvalue_sgi_z")
        ] = self._sgi_engine._pvalue_z
        
        data.iat[
            0, data.columns.get_loc("zcal")
        ] = self._sgi_engine._z_cal

        data.iat[
            0, data.columns.get_loc("pvalue_sgi_mc")
        ] = self._sgi_engine._pvalue_mc
        

        data.iat[
            0, data.columns.get_loc("Time (seconds) p-value SGI")
        ] = self._sgi_engine._time_pvalue

        data.fillna("", inplace=True)
        self._report_pvalue= data



    



    @property
    def report(self):
        self._build_report()

        return self._report

    
    @property
    def report_pvalue(self):
        if self._sgi_engine._pvalue_z is not None and self._sgi_engine._report_pvalue is None:
            self._build_report_pvalue()

        return self._report_pvalue