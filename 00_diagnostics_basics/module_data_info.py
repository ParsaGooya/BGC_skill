'''
class to store data information (e.g., obs, assimilation runs, simulations)

example:

class data_module(object):
  def __init__(self, loc):
    self.loc = loc
    
  def PrintLoc(self):
    print(self.loc)

a = data_module(loc)
a.PrintLoc()
'''


class data_module_chlos(object):

    def __init__(self, loc):

        self.loc = loc
        
        self.data_dict = { "ESACCI" : {"var"  : "chlos",
                                      "y0"    : 1998, 
                                      "y1"    : 2024,  
                                      "dir"   : f"{loc}" + "/ESACCI/ESACCI-1M_MONTHLY_1x1_199801-202403.nc",
                                      "file"  : "ESACCI-1M_MONTHLY_1x1",
                                      "color" : "black",
                                      "notes" : "ESACCI"
                                      },
                          "Shujie23" : {"var"  : "chlos",
                                      "y0"    : 1998, 
                                      "y1"    : 2020,  
                                      "dir"   : f"{loc}" + "/Shujie23/Shujie23-1M_MONTHLY_1x1_199709-202012.nc",
                                      "file"  : "Shujie23-1M_MONTHLY_1x1",
                                      "color" : "black",
                                      "notes" : "Shujie et al. (2023)"
                                      },
                           "asm"    : {"var"   : "chlos",
                                      "y0"    : 1958, 
                                      "y1"    : 2023,  
                                      "dir"   : f"{loc}" + "/chlos_Omon_ensmebles_*_1x1_LE.nc",
                                      "file"  : "A1-assim-re",
                                      "color" : "tab:green",
                                      "notes" : "None"
                                      },

                            "asm_CanOE"    : {"var"   : "chlos",
                                      "y0"    : 1958, 
                                      "y1"    : 2023,  
                                      "dir"   : f"{loc}" + "/chlos_Omon_ensmebles_*_1x1_LE.nc",
                                      "file"  : "A1-assim-re",
                                      "color" : "tab:cyan",
                                      "notes" : "None"
                                      },
                           "sim"   : {"var"   : "chlos",
                                      "y0"    : 1982, 
                                      "y1"    : 2024,  
                                      "dir"   : f"{loc}" + "/chlos_Omon_ensmebles_*_1x1_LE.nc",
                                      "file"  : "U1-simulation",
                                      "color" : "tab:red",
                                      "notes" : "None"
                                     },
                          "sim_CanOE"   : {"var"   : "chlos",
                                      "y0"    : 1950, 
                                      "y1"    : 2024,  
                                      "dir"   : f"{loc}" + "/chlos_Omon_ensmebles_*_1x1_LE.nc",
                                      "file"  : "U1-simulation",
                                      "color" : "tab:purple",
                                      "notes" : "None"
                                     },
                           "hnd"   : {"var"   : "fgco2",
                                      "y0"    : 1991, 
                                      "y1"    : 2024,  
                                      "dir"   : f"{loc}",
                                      "file"  : "fgco2_Omon",
                                      "color" : "tab:blue",
                                      "notes" : "None"
                                     },
                          #  "hnd_badj" : {"var"   : "fgco2",
                          #                "y0"    : 1987, 
                          #                "y1"    : 2021,  
                          #                "dir"   : f"{loc}/fgco2_ems/SOMFFN/results/Bias_Adjusted",
                          #                "file"  : "bias_adjusted_1987-2021",
                          #                "color" : "tab:cyan",
                          #                "notes" : "Bias adjusted out-of-sample using code bias_and_trend_correction.ipynb"
                          #               },
                    }
        
    def PrintLoc(self):
        print(self.loc)
    


class data_module_chl(object):

    def __init__(self, loc):

        self.loc = loc
        
        self.data_dict = { #"ESACCI" : {"var"  : "chl",
        #                               "y0"    : 1998, 
        #                               "y1"    : 2024,  
        #                               "dir"   : f"{loc}" + "/ESACCI/ESACCI-1M_MONTHLY_1x1_199801-202403.nc",
        #                               "file"  : "ESACCI-1M_MONTHLY_1x1",
        #                               "color" : "black",
        #                               "notes" : "ESACCI"
        #                               },
        #                   "Shujie23" : {"var"  : "chlos",
        #                               "y0"    : 1998, 
        #                               "y1"    : 2020,  
        #                               "dir"   : f"{loc}" + "/Shujie23/Shujie23-1M_MONTHLY_1x1_199709-202012.nc",
        #                               "file"  : "Shujie23-1M_MONTHLY_1x1",
        #                               "color" : "black",
        #                               "notes" : "Shujie et al. (2023)"
        #                               },
                           "asm"    : {"var"   : "chl",
                                      "y0"    : 1958, 
                                      "y1"    : 2023,  
                                      "dir"   : f"{loc}" + "/chl_Omon_ensmebles_*_1x1_LE.nc",
                                      "file"  : "A1-assim-re",
                                      "color" : "tab:green",
                                      "notes" : "None"
                                      },

                            "asm_CanOE"    : {"var"   : "chl",
                                      "y0"    : 1958, 
                                      "y1"    : 2023,  
                                      "dir"   : f"{loc}" + "/chl_Omon_ensmebles_*_1x1_LE.nc",
                                      "file"  : "A1-assim-re",
                                      "color" : "tab:cyan",
                                      "notes" : "None"
                                      },
                           "sim"   : {"var"   : "chl",
                                      "y0"    : 1982, 
                                      "y1"    : 2024,  
                                      "dir"   : f"{loc}" + "/chl_Omon_ensmebles_*_1x1_LE.nc",
                                      "file"  : "U1-simulation",
                                      "color" : "tab:red",
                                      "notes" : "None"
                                     },
                          "sim_CanOE"   : {"var"   : "chl",
                                      "y0"    : 1950, 
                                      "y1"    : 2024,  
                                      "dir"   : f"{loc}" + "/chl_Omon_ensmebles_*_1x1_LE.nc",
                                      "file"  : "U1-simulation",
                                      "color" : "tab:purple",
                                      "notes" : "None"
                                     },
                          #  "hnd"   : {"var"   : "fgco2",
                          #             "y0"    : 1991, 
                          #             "y1"    : 2024,  
                          #             "dir"   : f"{loc}",
                          #             "file"  : "fgco2_Omon",
                          #             "color" : "tab:blue",
                          #             "notes" : "None"
                          #            },
                          #  "hnd_badj" : {"var"   : "fgco2",
                          #                "y0"    : 1987, 
                          #                "y1"    : 2021,  
                          #                "dir"   : f"{loc}/fgco2_ems/SOMFFN/results/Bias_Adjusted",
                          #                "file"  : "bias_adjusted_1987-2021",
                          #                "color" : "tab:cyan",
                          #                "notes" : "Bias adjusted out-of-sample using code bias_and_trend_correction.ipynb"
                          #               },
                    }
        
    def PrintLoc(self):
        print(self.loc)


class data_module_pp(object):

    def __init__(self, loc):

        self.loc = loc
        
        self.data_dict = { "ESACCI" : {"var"  : "npp",
                                      "y0"    : 2002, 
                                      "y1"    : 2006,  
                                      "dir"   : f"{loc}" + "/cafe_global_1deg_2002_2006.nc",
                                      "file"  : "cafe_global_1deg_2002_2006",
                                      "color" : "black",
                                      "notes" : "ref?"
                                      },

                           "asm"    : {"var"   : "intpp",
                                      "y0"    : 1958, 
                                      "y1"    : 2020,  
                                      "dir"   : f"{loc}" + "/intpp_Omon_ensmebles_*_1x1_LE.nc",
                                      "file"  : "A1-assim-re",
                                      "color" : "tab:green",
                                      "notes" : "None"
                                      },
                           "sim"   : {"var"   : "intpp",
                                      "y0"    : 1980, 
                                      "y1"    : 2014,  
                                      "dir"   : f"{loc}" + "/intpp_Omon_ensmebles_*_1x1_LE.nc",
                                      "file"  : "U1-simulation",
                                      "color" : "tab:red",
                                      "notes" : "None"
                                     },
                           "hnd"   : {"var"   : "intpp",
                                      "y0"    : 1991, 
                                      "y1"    : 2020,  
                                      "dir"   : f"{loc}",
                                      "file"  : "int_Omon",
                                      "color" : "tab:blue",
                                      "notes" : "None"
                                     },

                    }
        
    def PrintLoc(self):
        print(self.loc)   



class data_module_no3(object):

    def __init__(self, loc):

        self.loc = loc
        
        self.data_dict = { "WOA23" : {"var"  : "no3",
                                      "y0"    : 1965, 
                                      "y1"    : 2022,  
                                      "dir"   : f"{loc}" + "/WOA23_clim_1x1_1965-2022.nc",
                                      "file"  : "WOA23_clim_1x1",
                                      "color" : "black",
                                      "notes" : "WOA23"
                                      },
                           "asm"    : {"var"   : "no3",
                                      "y0"    : 1958, 
                                      "y1"    : 2023,  
                                      "dir"   : f"{loc}" + "/no3_Omon_ensmebles_*_1x1_LE.nc",
                                      "file"  : "A1-assim-re",
                                      "color" : "tab:green",
                                      "notes" : "None"
                                      },

                            "asm_CanOE"    : {"var"   : "no3",
                                      "y0"    : 1958, 
                                      "y1"    : 2023,  
                                      "dir"   : f"{loc}" + "/no3_Omon_ensmebles_*_1x1_LE.nc",
                                      "file"  : "A1-assim-re",
                                      "color" : "tab:cyan",
                                      "notes" : "None"
                                      },
                           "sim"   : {"var"   : "no3",
                                      "y0"    : 1950, 
                                      "y1"    : 2025,  
                                      "dir"   : f"{loc}" + "/no3_Omon_ensmebles_*_1x1_LE.nc",
                                      "file"  : "U1-simulation",
                                      "color" : "tab:red",
                                      "notes" : "None"
                                     },
                          "sim_CanOE"   : {"var"   : "no3",
                                      "y0"    : 1950, 
                                      "y1"    : 2025,  
                                      "dir"   : f"{loc}" + "/no3_Omon_ensmebles_*_1x1_LE.nc",
                                      "file"  : "U1-simulation",
                                      "color" : "tab:purple",
                                      "notes" : "None"
                                     },
                          #  "hnd"   : {"var"   : "fgco2",
                          #             "y0"    : 1991, 
                          #             "y1"    : 2024,  
                          #             "dir"   : f"{loc}",
                          #             "file"  : "fgco2_Omon",
                          #             "color" : "tab:blue",
                          #             "notes" : "None"
                          #            },

                    }
        
    def PrintLoc(self):
        print(self.loc)
        

class data_module_wo(object):

    def __init__(self, loc):

        self.loc = loc
        
        self.data_dict = { 
                          "obs" : {"var"  : "wo",
                                      "y0"    : 1980, 
                                      "y1"    : 2023,  
                                      "dir"   : f"{loc}" + "/wo_Omon_obs_198001_202312_1x1.nc",
                                      "file"  : "soda_1x1",
                                      "color" : "black",
                                      "notes" : "soda_3.15.2"
                                      },
                           "asm"    : {"var"   : "wo",
                                      "y0"    : 1958, 
                                      "y1"    : 2023,  
                                      "dir"   : f"{loc}" + "/wo_Omon_ensmebles_*_1x1_LE.nc",
                                      "file"  : "A1-assim-re",
                                      "color" : "tab:green",
                                      "notes" : "None"
                                      },


                           "sim"   : {"var"   : "wo",
                                      "y0"    : 1950, 
                                      "y1"    : 2025,  
                                      "dir"   : f"{loc}" + "/wo_Omon_ensmebles_*_1x1_LE.nc",
                                      "file"  : "U1-simulation",
                                      "color" : "tab:red",
                                      "notes" : "None"
                                     },

                          #  "hnd"   : {"var"   : "fgco2",
                          #             "y0"    : 1991, 
                          #             "y1"    : 2024,  
                          #             "dir"   : f"{loc}",
                          #             "file"  : "fgco2_Omon",
                          #             "color" : "tab:blue",
                          #             "notes" : "None"
                          #            },

                    }
        
    def PrintLoc(self):
        print(self.loc)


class data_module_thetao(object):

    def __init__(self, loc):

        self.loc = loc
        
        self.data_dict = { 
                          "obs" : {"var"  : "thetao",
                                      "y0"    : 1980, 
                                      "y1"    : 2023,  
                                      "dir"   : f"{loc}" + "/thetao_Omon_obs_198001_202312_1x1.nc",
                                      "file"  : "soda_1x1",
                                      "color" : "black",
                                      "notes" : "soda_3.15.2"
                                      },
                           "asm"    : {"var"   : "thetao",
                                      "y0"    : 1958, 
                                      "y1"    : 2023,  
                                      "dir"   : f"{loc}" + "/thetao_Omon_ensmebles_*_1x1_LE.nc",
                                      "file"  : "A1-assim-re",
                                      "color" : "tab:green",
                                      "notes" : "None"
                                      },


                           "sim"   : {"var"   : "thetao",
                                      "y0"    : 1950, 
                                      "y1"    : 2025,  
                                      "dir"   : f"{loc}" + "/thetao_Omon_ensmebles_*_1x1_LE.nc",
                                      "file"  : "U1-simulation",
                                      "color" : "tab:red",
                                      "notes" : "None"
                                     },

                          #  "hnd"   : {"var"   : "fgco2",
                          #             "y0"    : 1991, 
                          #             "y1"    : 2024,  
                          #             "dir"   : f"{loc}",
                          #             "file"  : "fgco2_Omon",
                          #             "color" : "tab:blue",
                          #             "notes" : "None"
                          #            },

                    }
        
    def PrintLoc(self):
        print(self.loc)




class data_module_uo(object):

    def __init__(self, loc):

        self.loc = loc
        
        self.data_dict = { 
                          "obs" : {"var"  : "uo",
                                      "y0"    : 1980, 
                                      "y1"    : 2023,  
                                      "dir"   : f"{loc}" + "/uo_Omon_obs_198001_202312_1x1.nc",
                                      "file"  : "soda_1x1",
                                      "color" : "black",
                                      "notes" : "soda_3.15.2"
                                      },
                           "asm"    : {"var"   : "uo",
                                      "y0"    : 1958, 
                                      "y1"    : 2023,  
                                      "dir"   : f"{loc}" + "/uo_Omon_ensmebles_*_1x1_LE.nc",
                                      "file"  : "A1-assim-re",
                                      "color" : "tab:green",
                                      "notes" : "None"
                                      },


                           "sim"   : {"var"   : "uo",
                                      "y0"    : 1950, 
                                      "y1"    : 2025,  
                                      "dir"   : f"{loc}" + "/uo_Omon_ensmebles_*_1x1_LE.nc",
                                      "file"  : "U1-simulation",
                                      "color" : "tab:red",
                                      "notes" : "None"
                                     },

                          #  "hnd"   : {"var"   : "fgco2",
                          #             "y0"    : 1991, 
                          #             "y1"    : 2024,  
                          #             "dir"   : f"{loc}",
                          #             "file"  : "fgco2_Omon",
                          #             "color" : "tab:blue",
                          #             "notes" : "None"
                          #            },

                    }
        
    def PrintLoc(self):
        print(self.loc)



class data_module_vo(object):

    def __init__(self, loc):

        self.loc = loc
        
        self.data_dict = { 
                          "obs" : {"var"  : "vo",
                                      "y0"    : 1980, 
                                      "y1"    : 2023,  
                                      "dir"   : f"{loc}" + "/vo_Omon_obs_198001_202312_1x1.nc",
                                      "file"  : "soda_1x1",
                                      "color" : "black",
                                      "notes" : "soda_3.15.2"
                                      },
                           "asm"    : {"var"   : "vo",
                                      "y0"    : 1958, 
                                      "y1"    : 2023,  
                                      "dir"   : f"{loc}" + "/vo_Omon_ensmebles_*_1x1_LE.nc",
                                      "file"  : "A1-assim-re",
                                      "color" : "tab:green",
                                      "notes" : "None"
                                      },


                           "sim"   : {"var"   : "vo",
                                      "y0"    : 1950, 
                                      "y1"    : 2025,  
                                      "dir"   : f"{loc}" + "/vo_Omon_ensmebles_*_1x1_LE.nc",
                                      "file"  : "U1-simulation",
                                      "color" : "tab:red",
                                      "notes" : "None"
                                     },

                          #  "hnd"   : {"var"   : "fgco2",
                          #             "y0"    : 1991, 
                          #             "y1"    : 2024,  
                          #             "dir"   : f"{loc}",
                          #             "file"  : "fgco2_Omon",
                          #             "color" : "tab:blue",
                          #             "notes" : "None"
                          #            },

                    }
        
    def PrintLoc(self):
        print(self.loc)


class data_module_so(object):

    def __init__(self, loc):

        self.loc = loc
        
        self.data_dict = { 
                          "obs" : {"var"  : "so",
                                      "y0"    : 1980, 
                                      "y1"    : 2023,  
                                      "dir"   : f"{loc}" + "/so_Omon_obs_198001_202312_1x1.nc",
                                      "file"  : "soda_1x1",
                                      "color" : "black",
                                      "notes" : "soda_3.15.2"
                                      },
                           "asm"    : {"var"   : "so",
                                      "y0"    : 1958, 
                                      "y1"    : 2023,  
                                      "dir"   : f"{loc}" + "/so_Omon_ensmebles_*_1x1_LE.nc",
                                      "file"  : "A1-assim-re",
                                      "color" : "tab:green",
                                      "notes" : "None"
                                      },


                           "sim"   : {"var"   : "so",
                                      "y0"    : 1950, 
                                      "y1"    : 2025,  
                                      "dir"   : f"{loc}" + "/so_Omon_ensmebles_*_1x1_LE.nc",
                                      "file"  : "U1-simulation",
                                      "color" : "tab:red",
                                      "notes" : "None"
                                     },

                          #  "hnd"   : {"var"   : "fgco2",
                          #             "y0"    : 1991, 
                          #             "y1"    : 2024,  
                          #             "dir"   : f"{loc}",
                          #             "file"  : "fgco2_Omon",
                          #             "color" : "tab:blue",
                          #             "notes" : "None"
                          #            },

                    }
        
    def PrintLoc(self):
        print(self.loc)


class data_module_agessc(object):

    def __init__(self, loc):

        self.loc = loc
        
        self.data_dict = { 
                          # "obs" : {"var"  : "so",
                          #             "y0"    : 1980, 
                          #             "y1"    : 2023,  
                          #             "dir"   : f"{loc}" + "/agessc_Omon_obs_198001_202312_1x1.nc",
                          #             "file"  : "soda_1x1",
                          #             "color" : "black",
                          #             "notes" : "soda_3.15.2"
                          #             },
                           "asm"    : {"var"   : "agessc",
                                      "y0"    : 1958, 
                                      "y1"    : 2023,  
                                      "dir"   : f"{loc}" + "/agessc_Omon_ensmebles_*_1x1_LE.nc",
                                      "file"  : "A1-assim-re",
                                      "color" : "tab:green",
                                      "notes" : "None"
                                      },


                           "sim"   : {"var"   : "agessc",
                                      "y0"    : 1950, 
                                      "y1"    : 2025,  
                                      "dir"   : f"{loc}" + "/agessc_Omon_ensmebles_*_1x1_LE.nc",
                                      "file"  : "U1-simulation",
                                      "color" : "tab:red",
                                      "notes" : "None"
                                     },

                          #  "hnd"   : {"var"   : "fgco2",
                          #             "y0"    : 1991, 
                          #             "y1"    : 2024,  
                          #             "dir"   : f"{loc}",
                          #             "file"  : "fgco2_Omon",
                          #             "color" : "tab:blue",
                          #             "notes" : "None"
                          #            },

                    }
        
    def PrintLoc(self):
        print(self.loc)

class data_module_uas(object):

    def __init__(self, loc):

        self.loc = loc
        
        self.data_dict = { 
                          "obs" : {"var"  : "uas",
                                      "y0"    : 1980, 
                                      "y1"    : 2023,  
                                      "dir"   : f"{loc}" + "/uas_Amon_obs_*_1x1.nc",
                                      "file"  : "ERA5_1x1",
                                      "color" : "black",
                                      "notes" : "soda_3.15.2"
                                      },
                           "asm"    : {"var"   : "uas",
                                      "y0"    : 1958, 
                                      "y1"    : 2023,  
                                      "dir"   : f"{loc}" + "/uas_Amon_ensmebles_*_1x1_LE.nc",
                                      "file"  : "A1-assim-re",
                                      "color" : "tab:green",
                                      "notes" : "None"
                                      },


                           "sim"   : {"var"   : "uas",
                                      "y0"    : 1950, 
                                      "y1"    : 2025,  
                                      "dir"   : f"{loc}" + "/uas_Amon_ensmebles_*_1x1_LE.nc",
                                      "file"  : "U1-simulation",
                                      "color" : "tab:red",
                                      "notes" : "None"
                                     },

                          #  "hnd"   : {"var"   : "fgco2",
                          #             "y0"    : 1991, 
                          #             "y1"    : 2024,  
                          #             "dir"   : f"{loc}",
                          #             "file"  : "fgco2_Omon",
                          #             "color" : "tab:blue",
                          #             "notes" : "None"
                          #            },

                    }
        
    def PrintLoc(self):
        print(self.loc)

class data_module_vas(object):

    def __init__(self, loc):

        self.loc = loc
        
        self.data_dict = { 
                          "obs" : {"var"  : "vas",
                                      "y0"    : 1980, 
                                      "y1"    : 2023,  
                                      "dir"   : f"{loc}" + "/vas_Amon_obs_*_1x1.nc",
                                      "file"  : "ERA5_1x1",
                                      "color" : "black",
                                      "notes" : "soda_3.15.2"
                                      },
                           "asm"    : {"var"   : "vas",
                                      "y0"    : 1958, 
                                      "y1"    : 2023,  
                                      "dir"   : f"{loc}" + "/vas_Amon_ensmebles_*_1x1_LE.nc",
                                      "file"  : "A1-assim-re",
                                      "color" : "tab:green",
                                      "notes" : "None"
                                      },


                           "sim"   : {"var"   : "vas",
                                      "y0"    : 1950, 
                                      "y1"    : 2025,  
                                      "dir"   : f"{loc}" + "/vas_Amon_ensmebles_*_1x1_LE.nc",
                                      "file"  : "U1-simulation",
                                      "color" : "tab:red",
                                      "notes" : "None"
                                     },

                          #  "hnd"   : {"var"   : "fgco2",
                          #             "y0"    : 1991, 
                          #             "y1"    : 2024,  
                          #             "dir"   : f"{loc}",
                          #             "file"  : "fgco2_Omon",
                          #             "color" : "tab:blue",
                          #             "notes" : "None"
                          #            },

                    }
        
    def PrintLoc(self):
        print(self.loc)


class biomes_module(object):
    
    def __init__(self, loc):
        
        self.loc = loc
        
        self.data_dict = {
                            "biomes" : {
                                        "y0"         : 1998, 
                                        "y1"         : 2010,  
                                        "dir"        : f"{loc}",
                                        "file"       : "Time_Varying_Biomes",
                                        "area_units" : "1e6 km2",
                                        "area_ocean" : 361.0,
                                        "notes" : "Global Ocean Biomes as per Fay and McKinley (2014)), downloaded from: doi:10.1594/PANGAEA.828650"
                                      },            
                          }
        
        self.biomes_dict = {
                 1 : { 'label'     : 'NP_ICE',
                       'mean_area' : 4.5852,
                       'core_area' : 3.9112,
                     },
                 2 : { 'label'     : 'NP_SPSS',
                       'mean_area' : 12.838,
                       'core_area' : 8.7063,
                     },
                 3 : { 'label'     : 'NP_STSS',
                       'mean_area' : 6.8257,
                       'core_area' : 4.3398,
                     },
                 4 : { 'label'     : 'NP_STPS',
                       'mean_area' : 41.048,
                       'core_area' : 31.305,
                     },
                 5 : { 'label'     : 'P_EQU_W', #'Pac_EQU_W',
                       'mean_area' : 11.593,
                       'core_area' : 4.9040,
                     },
                 6 : { 'label'     : 'P_EQU_E',#'Pac_EQU_E',
                       'mean_area' : 14.890,
                       'core_area' : 7.7117,
                     },
                 7 : { 'label'     : 'SP_STPS',
                       'mean_area' : 52.705,
                       'core_area' : 43.067,
                     },
                 8 : { 'label'     : 'NA_ICE',
                       'mean_area' : 5.4750,
                       'core_area' : 4.6123,
                     },
                 9 : { 'label'     : 'NA_SPSS',
                       'mean_area' : 10.062,
                       'core_area' : 7.8048,
                     },
                 10 : { 'label'    : 'NA_STSS',
                       'mean_area' : 5.9744,
                       'core_area' : 4.6178,
                     },
                 11 : { 'label'    : 'NA_STPS',
                       'mean_area' : 17.464,
                       'core_area' : 13.746,
                     },
                 12 : { 'label'    : 'A_EQU', #'Atl_EQU',
                       'mean_area' : 7.4147,
                       'core_area' : 3.1972,
                     },
                 13 : { 'label'    : 'SA_STPS',
                       'mean_area' : 18.055,
                       'core_area' : 16.056,
                     },
                 14 : { 'label'    : 'IND_STPS',
                       'mean_area' : 35.936,
                       'core_area' : 32.383,
                     },
                 15 : { 'label'    : 'SO_STSS',
                       'mean_area' : 29.692,
                       'core_area' : 23.670,
                     },
                 16 : { 'label'    : 'SO_SPSS',
                       'mean_area' : 30.628,
                       'core_area' : 25.875,
                     },
                 17 : { 'label'    : 'SO_ICE',
                       'mean_area' : 18.678,
                       'core_area' : 16.169,
                     }
               }    
        
        self.dict_biomes_plot = {
                 1 : { 'label'        : 'NP_ICE',
                       'ts_clim_mean_ymin'    : -.05,
                       'ts_clim_mean_ymax'    : 1.,
                       'ts_clim_var_ymin'     : 0,
                       'ts_clim_var_ymax'     : .3,
                       'ts_avg_ymin'        : -1e12,
                       'ts_avg_ymax'        : 4e12,
                       'ts_avg_rmse_ymin'   : 0.,
                       'ts_avg_rmse_ymax'   : 10.,
                       'ts_avg_corr_ymin'   : 0., 
                       'ts_avg_corr_ymax'   : .7, 
                       'ts_rmse_avg_ymin'     : 0.4,
                       'ts_rmse_avg_ymax'     : 1.,
                       'ts_corr_avg_ymin'     : 0.1, 
                       'ts_corr_avg_ymax'     : 0.5,
                     },
                 2 : { 'label'     : 'NP_SPSS',
                       'ts_clim_mean_ymin'    : -.05,
                       'ts_clim_mean_ymax'    : 1.,
                       'ts_clim_var_ymin'     : 0,
                       'ts_clim_var_ymax'     : .3,
                       'ts_avg_ymin'        : -1e12,
                       'ts_avg_ymax'        : 4e12,
                       'ts_avg_rmse_ymin'   : 0.,
                       'ts_avg_rmse_ymax'   : 10.,
                       'ts_avg_corr_ymin'   : 0., 
                       'ts_avg_corr_ymax'   : .7, 
                       'ts_rmse_avg_ymin'     : 0.3,
                       'ts_rmse_avg_ymax'     : 1.5,
                       'ts_corr_avg_ymin'     : 0.1, 
                       'ts_corr_avg_ymax'     : 0.7,

                     },
                 3 : { 'label'     : 'NP_STSS',
                       'ts_clim_mean_ymin'    : -.05,
                       'ts_clim_mean_ymax'    : 1.,
                       'ts_clim_var_ymin'     : 0,
                       'ts_clim_var_ymax'     : .3,
                       'ts_avg_ymin'        : -1e12,
                       'ts_avg_ymax'        : 4e12,
                       'ts_avg_rmse_ymin'   : 0.,
                       'ts_avg_rmse_ymax'   : 10.,
                       'ts_avg_corr_ymin'   : 0., 
                       'ts_avg_corr_ymax'   : .7, 
                       'ts_rmse_avg_ymin'     : 0.,
                       'ts_rmse_avg_ymax'     : 2.5,
                       'ts_corr_avg_ymin'     : -0.3, 
                       'ts_corr_avg_ymax'     : 0.5,
                     },
                 4 : { 'label'     : 'NP_STPS',
                       'ts_clim_mean_ymin'    : -.05,
                       'ts_clim_mean_ymax'    : 1.,
                       'ts_clim_var_ymin'     : 0,
                       'ts_clim_var_ymax'     : .3,
                       'ts_avg_ymin'        : -1e12,
                       'ts_avg_ymax'        : 4e12,
                       'ts_avg_rmse_ymin'   : 0.,
                       'ts_avg_rmse_ymax'   : 10.,
                       'ts_avg_corr_ymin'   : 0., 
                       'ts_avg_corr_ymax'   : .7, 
                       'ts_rmse_avg_ymin'     : 0.1,
                       'ts_rmse_avg_ymax'     : .5,
                       'ts_corr_avg_ymin'     : -.1, 
                       'ts_corr_avg_ymax'     : 0.4,
                     },
                 5 : { 'label'     : 'Pac_EQU_W',
                       'ts_clim_mean_ymin'    : -.05,
                       'ts_clim_mean_ymax'    : 1.,
                       'ts_clim_var_ymin'     : 0,
                       'ts_clim_var_ymax'     : .3,
                       'ts_avg_ymin'        : -1e12,
                       'ts_avg_ymax'        : 4e12,
                       'ts_avg_rmse_ymin'   : 0.,
                       'ts_avg_rmse_ymax'   : 10.,
                       'ts_avg_corr_ymin'   : 0., 
                       'ts_avg_corr_ymax'   : .7, 
                       'ts_rmse_avg_ymin'     : 0.,
                       'ts_rmse_avg_ymax'     : .7,
                       'ts_corr_avg_ymin'     : 0., 
                       'ts_corr_avg_ymax'     : 0.5,
                     },
                 6 : { 'label'     : 'Pac_EQU_E',
                       'ts_clim_mean_ymin'    : 1.,
                       'ts_clim_mean_ymax'    : 3.5,
                       'ts_clim_var_ymin'     : 0.,
                       'ts_clim_var_ymax'     : 0.8,
                       'ts_avg_ymin'        : 0.9e13,
                       'ts_avg_ymax'        : 2.5e13,
                       'ts_avg_rmse_ymin'   : 0.1, 
                       'ts_avg_rmse_ymax'   : 0.4,
                       'ts_avg_corr_ymin'   : -.4, 
                       'ts_avg_corr_ymax'   : .8, 
                       'ts_rmse_avg_ymin'     : .3,
                       'ts_rmse_avg_ymax'     : 1.1,
                       'ts_corr_avg_ymin'     : -0.2, 
                       'ts_corr_avg_ymax'     : 0.6, 
                     },
                 7 : { 'label'     : 'SP_STPS',
                       'ts_clim_mean_ymin'    : -.05,
                       'ts_clim_mean_ymax'    : 1.,
                       'ts_clim_var_ymin'     : 0,
                       'ts_clim_var_ymax'     : .3,
                       'ts_avg_ymin'        : -1e12,
                       'ts_avg_ymax'        : 4e12,
                       'ts_avg_rmse_ymin'   : 0.,
                       'ts_avg_rmse_ymax'   : 10.,
                       'ts_avg_corr_ymin'   : 0., 
                       'ts_avg_corr_ymax'   : .7, 
                       'ts_rmse_avg_ymin'     : 0.,
                       'ts_rmse_avg_ymax'     : 1.,
                       'ts_corr_avg_ymin'     : 0.2, 
                       'ts_corr_avg_ymax'     : 0.6,
                     },
                 8 : { 'label'     : 'NA_ICE',
                       'ts_clim_mean_ymin'    : -.05,
                       'ts_clim_mean_ymax'    : 1.,
                       'ts_clim_var_ymin'     : 0,
                       'ts_clim_var_ymax'     : .3,
                       'ts_avg_ymin'        : -1e12,
                       'ts_avg_ymax'        : 4e12,
                       'ts_avg_rmse_ymin'   : 0.,
                       'ts_avg_rmse_ymax'   : 10.,
                       'ts_avg_corr_ymin'   : 0., 
                       'ts_avg_corr_ymax'   : .7, 
                       'ts_rmse_avg_ymin'     : 0.2,
                       'ts_rmse_avg_ymax'     : 1.6,
                       'ts_corr_avg_ymin'     : 0.2, 
                       'ts_corr_avg_ymax'     : 0.5,
                     },
                 9 : { 'label'     : 'NA_SPSS',
                       'ts_clim_mean_ymin'    : -.05,
                       'ts_clim_mean_ymax'    : 1.,
                       'ts_clim_var_ymin'     : 0,
                       'ts_clim_var_ymax'     : .3,
                       'ts_avg_ymin'        : -1e12,
                       'ts_avg_ymax'        : 4e12,
                       'ts_avg_rmse_ymin'   : 0.,
                       'ts_avg_rmse_ymax'   : 10.,
                       'ts_avg_corr_ymin'   : 0., 
                       'ts_avg_corr_ymax'   : .7, 
                       'ts_rmse_avg_ymin'     : 0.4,
                       'ts_rmse_avg_ymax'     : 1.8,
                       'ts_corr_avg_ymin'     : 0., 
                       'ts_corr_avg_ymax'     : 0.8,
                     },
                 10 : { 'label'     : 'NA_STSS',
                       'ts_clim_mean_ymin'    : -.05,
                       'ts_clim_mean_ymax'    : 1.,
                       'ts_clim_var_ymin'     : 0,
                       'ts_clim_var_ymax'     : .3,
                       'ts_avg_ymin'        : -1e12,
                       'ts_avg_ymax'        : 4e12,
                       'ts_avg_rmse_ymin'   : 0.,
                       'ts_avg_rmse_ymax'   : 10.,
                       'ts_avg_corr_ymin'   : 0., 
                       'ts_avg_corr_ymax'   : .7, 
                       'ts_rmse_avg_ymin'     : 0.1,
                       'ts_rmse_avg_ymax'     : 1.6,
                       'ts_corr_avg_ymin'     : 0.2, 
                       'ts_corr_avg_ymax'     : 0.8,
                     },
                 11 : { 'label'    : 'NA_STPS',
                       'ts_clim_mean_ymin'    : -.05,
                       'ts_clim_mean_ymax'    : 1.,
                       'ts_clim_var_ymin'     : 0,
                       'ts_clim_var_ymax'     : .3,
                       'ts_avg_ymin'        : -1e12,
                       'ts_avg_ymax'        : 4e12,
                       'ts_avg_rmse_ymin'   : 0.,
                       'ts_avg_rmse_ymax'   : 10.,
                       'ts_avg_corr_ymin'   : 0., 
                       'ts_avg_corr_ymax'   : .7, 
                       'ts_rmse_avg_ymin'     : 0.1,
                       'ts_rmse_avg_ymax'     : 0.4,
                       'ts_corr_avg_ymin'     : 0., 
                       'ts_corr_avg_ymax'     : 0.5,
                     },
                 12 : { 'label'    : 'Atl_EQU',
                       'ts_clim_mean_ymin'    : -.05,
                       'ts_clim_mean_ymax'    : 1.,
                       'ts_clim_var_ymin'     : 0,
                       'ts_clim_var_ymax'     : .3,
                       'ts_avg_ymin'        : -1e12,
                       'ts_avg_ymax'        : 4e12,
                       'ts_avg_rmse_ymin'   : 0.,
                       'ts_avg_rmse_ymax'   : 10.,
                       'ts_avg_corr_ymin'   : 0., 
                       'ts_avg_corr_ymax'   : .7, 
                       'ts_rmse_avg_ymin'     : 0.,
                       'ts_rmse_avg_ymax'     : .4,
                       'ts_corr_avg_ymin'     : 0., 
                       'ts_corr_avg_ymax'     : 0.7,
                     },
                 13 : { 'label'    : 'SA_STPS',
                       'ts_clim_mean_ymin'    : -.05,
                       'ts_clim_mean_ymax'    : 1.,
                       'ts_clim_var_ymin'     : 0,
                       'ts_clim_var_ymax'     : .3,
                       'ts_avg_ymin'        : -1e12,
                       'ts_avg_ymax'        : 4e12,
                       'ts_avg_rmse_ymin'   : 0.,
                       'ts_avg_rmse_ymax'   : 10.,
                       'ts_avg_corr_ymin'   : 0., 
                       'ts_avg_corr_ymax'   : .7, 
                       'ts_rmse_avg_ymin'     : 0.,
                       'ts_rmse_avg_ymax'     : 1.2,
                       'ts_corr_avg_ymin'     : 0., 
                       'ts_corr_avg_ymax'     : 0.6,
                     },
                 14 : { 'label'    : 'IND_STPS',
                       'ts_clim_mean_ymin'    : -.05,
                       'ts_clim_mean_ymax'    : 1.,
                       'ts_clim_var_ymin'     : 0,
                       'ts_clim_var_ymax'     : .3,
                       'ts_avg_ymin'        : -1e12,
                       'ts_avg_ymax'        : 4e12,
                       'ts_avg_rmse_ymin'   : 0.,
                       'ts_avg_rmse_ymax'   : 10.,
                       'ts_avg_corr_ymin'   : 0., 
                       'ts_avg_corr_ymax'   : .7, 
                       'ts_rmse_avg_ymin'     : 0.1,
                       'ts_rmse_avg_ymax'     : .6,
                       'ts_corr_avg_ymin'     : 0.1, 
                       'ts_corr_avg_ymax'     : 0.6,
                     },
                 15 : { 'label'     : 'SO_STSS',
                       'ts_clim_mean_ymin'    : -2.5,
                       'ts_clim_mean_ymax'    : 1.1,
                       'ts_clim_var_ymin'     : 0,
                       'ts_clim_var_ymax'     : 1.4,
                       'ts_avg_ymin'        : -.7e14, 
                       'ts_avg_ymax'        : 0.,
                       'ts_avg_rmse_ymin'   : 0.,
                       'ts_avg_rmse_ymax'   : 1.,
                       'ts_avg_corr_ymin'   : 0.4, 
                       'ts_avg_corr_ymax'   : 1., 
                       'ts_rmse_avg_ymin'     : 0.,
                       'ts_rmse_avg_ymax'     : 2., 
                       'ts_corr_avg_ymin'     : 0.2, 
                       'ts_corr_avg_ymax'     : 0.8,
                     },
                 16 : { 'label'    : 'SO_SPSS',
                       'ts_clim_mean_ymin'    : -1.5,
                       'ts_clim_mean_ymax'    : 2.5,
                       'ts_clim_var_ymin'     : 0,
                       'ts_clim_var_ymax'     : .6,
                       'ts_avg_ymin'        : -3e13,
                       'ts_avg_ymax'        : 2e13,
                       'ts_avg_rmse_ymin'   : 8,
                       'ts_avg_rmse_ymax'   : 33,
                       'ts_avg_corr_ymin'   : -.3, 
                       'ts_avg_corr_ymax'   : 1., 
                       'ts_rmse_avg_ymin'     : .2,
                       'ts_rmse_avg_ymax'     : 1.2,
                       'ts_corr_avg_ymin'     : 0., 
                       'ts_corr_avg_ymax'     : .8,
                     },
                 17 : { 'label'    : 'SO_ICE',
                       'ts_clim_mean_ymin'    : -.05,
                       'ts_clim_mean_ymax'    : 1.,
                       'ts_clim_var_ymin'     : 0,
                       'ts_clim_var_ymax'     : .3,
                       'ts_avg_ymin'        : -1e12,
                       'ts_avg_ymax'        : 4e12,
                       'ts_avg_rmse_ymin'   : 0.,
                       'ts_avg_rmse_ymax'   : 10.,
                       'ts_avg_corr_ymin'   : 0., 
                       'ts_avg_corr_ymax'   : .7, 
                       'ts_rmse_avg_ymin'     : 0.2,
                       'ts_rmse_avg_ymax'     : 0.7,
                       'ts_corr_avg_ymin'     : -.1, 
                       'ts_corr_avg_ymax'     : 0.5,
                     }
               }    
    
    def PrintLoc(self):
        print(self.loc)

        
        
