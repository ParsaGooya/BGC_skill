

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

      
# # class data_modules(object):

# #     def __init__(self, loc, var, assimilation_BGC_run_id = None, CanOE_assimilation_BGC_run_id = 1):

# #         if assimilation_BGC_run_id is not None:
# #             if all([loc.split('/')[-1] == 'CanESM5','assimilation' in loc]):                
# #                 loc = loc.split('assimilation')[0] + f'assimilation_bgc{assimilation_BGC_run_id}' + loc.split('assimilation')[1]

# #         if all(['CanESM5-CanOE' in loc.split('/')[-1] ,'assimilation' in loc]):        
# #                 if CanOE_assimilation_BGC_run_id is not None:        
# #                   loc = loc.split('CanESM5-CanOE')[0] + f'CanESM5-CanOE_{CanOE_assimilation_BGC_run_id}'  
# #                 else:
# #                   raise RuntimeError("Please specify CanOE_assimilation_BGC_run_id ...")
# #         self.loc = loc

# #         if 'obs' in self.loc:
               
# #           if var == 'chlos':
# #                   self.data_dict = { "ESACCI" : {"var"  : "chlos",
# #                                             "y0"    : 1998, 
# #                                             "y1"    : 2024,  
# #                                             "dir"   : f"{loc}" + "/ESACCI/ESACCI-1M_MONTHLY_1x1_199801-202403.nc",
# #                                             "file"  : "ESACCI-1M_MONTHLY_1x1",
# #                                             "color" : "black",
# #                                             "notes" : "ESACCI"
# #                                             },
# #                                 "Shujie23" : {"var"  : "chlos",
# #                                             "y0"    : 1998, 
# #                                             "y1"    : 2020,  
# #                                             "dir"   : f"{loc}" + "/Shujie23/Shujie23-1M_MONTHLY_1x1_199709-202012.nc",
# #                                             "file"  : "Shujie23-1M_MONTHLY_1x1",
# #                                             "color" : "black",
# #                                             "notes" : "Shujie et al. (2023)"
# #                                             }}
# #           elif var == 'intpp':
# #                       self.data_dict = { "ESACCI" : {"var"  : "npp",
# #                                         "y0"    : 2002, 
# #                                         "y1"    : 2006,  
# #                                         "dir"   : f"{loc}" + "/cafe_global_1deg_2002_2006.nc",
# #                                         "file"  : "cafe_global_1deg_2002_2006",
# #                                         "color" : "black",
# #                                         "notes" : "ref?"
# #                                         } }     
# #           elif var == 'no3':
# #                      self.data_dict = { "WOA23" : {"var"  : "no3",
# #                                       "y0"    : 1965, 
# #                                       "y1"    : 2022,  
# #                                       "dir"   : f"{loc}" + "/WOA23_clim_1x1_1965-2022.nc",
# #                                       "file"  : "WOA23_clim_1x1",
# #                                       "color" : "black",
# #                                       "notes" : "WOA23"
# #                                       },
# #                             "GLODAP" : {"var"  : var,
# #                                       "y0"    : 1980, 
# #                                       "y1"    : 2023,  
# #                                       "dir"   : f"{loc}" + f"/GLODAP_v2_2023_{var}_1984-2021.csv",
# #                                       "file"  : "GLODAP",
# #                                       "color" : "black",
# #                                       "notes" : "GLODAP"
# #                                       },
# #                             "GLODAP_climatology" : {"var"  : var,
# #                                       "y0"    : 1980, 
# #                                       "y1"    : 2016,  
# #                                       "dir"   : f"{loc}" + f"/GLODAPv2.2016b.{var}.nc",
# #                                       "file"  : "GLODAP",
# #                                       "color" : "black",
# #                                       "notes" : "GLODAP"
# #                                       }}
# #           elif var in ['po4', 'silicate', 'o2']:
# #                      self.data_dict = { 
# #                             "GLODAP" : {"var"  : var,
# #                                       "y0"    : 1980, 
# #                                       "y1"    : 2023,  
# #                                       "dir"   : f"{loc}" + f"/GLODAP_v2_2023_{var}_1984-2021.csv",
# #                                       "file"  : "GLODAP",
# #                                       "color" : "black",
# #                                       "notes" : "GLODAP"
# #                                       },
# #                             "GLODAP_climatology" : {"var"  : var,
# #                                       "y0"    : 1980, 
# #                                       "y1"    : 2016,  
# #                                       "dir"   : f"{loc}" + f"/GLODAPv2.2016b.{var}.nc",
# #                                       "file"  : "GLODAP",
# #                                       "color" : "black",
# #                                       "notes" : "GLODAP"
# #                                       }}
# #           elif var == 'spco2':
             
# #                     self.data_dict = { "obs" : {"var"  : "spco2",
# #                                       "y0"    : 1982, 
# #                                       "y1"    : 2023,  
# #                                       "dir"   : f"{loc}" + "/spco2_ice_masked_MPI_SOMFFN_1x1_1982-2022.nc",
# #                                       "file"  : "MPI_SOM_FFN_1x1",
# #                                       "color" : "black",
# #                                       "notes" : "SOMFFNv2023"
# #                                       }}
# #           elif var == 'tos':

# #                 self.data_dict = { 
# #                           "obs" : {"var"  : "tos",
# #                                       "y0"    : 1980, 
# #                                       "y1"    : 2022,  
# #                                       "dir"   : f"{loc}" + "/ERSST_v5_tos_198001-202112_1x1.nc",
# #                                       "file"  : "ersst_v5",
# #                                       "color" : "black",
# #                                       }}
                
# #           elif var in  ['uas','vas','tauu','tauv']:
# #                      self.data_dict = { 
# #                           "obs" : {"var"  : var,
# #                                       "y0"    : 1980, 
# #                                       "y1"    : 2023,  
# #                                       "dir"   : f"{loc}" + f"/ERA5/{var}_Amon_obs_*_1x1.nc",
# #                                       "file"  : "ERA5_1x1",
# #                                       "color" : "black",
# #                                       "notes" : "soda_3.15.2"
# #                                       },}
# #           elif var in ['wo','uo','vo']:
# #                       self.data_dict = { 
# #                           "obs" : {"var"  : var,
# #                                       "y0"    : 1980, 
# #                                       "y1"    : 2023,  
# #                                       "dir"   : f"{loc}" + f"/SODA/{var}_Omon_obs_198001_202312_1x1.nc",
# #                                       "file"  : "soda_1x1",
# #                                       "color" : "black",
# #                                       "notes" : "soda_3.15.2"
# #                                       },}
# #           elif var in ['thetao','so']:
# #                       self.data_dict = { 
# #                           "SODA" : {"var"  : var,
# #                                       "y0"    : 1980, 
# #                                       "y1"    : 2023,  
# #                                       "dir"   : f"{loc}" + f"/SODA/{var}_Omon_obs_198001_202312_1x1.nc",
# #                                       "file"  : "soda_1x1",
# #                                       "color" : "black",
# #                                       "notes" : "soda_3.15.2"
# #                                       },
# #                             "GLODAP" : {"var"  : var,
# #                                       "y0"    : 1980, 
# #                                       "y1"    : 2023,  
# #                                       "dir"   : f"{loc}" + f"/GLODAP_v2_2023_{var}_1984-2021.csv",
# #                                       "file"  : "GLODAP",
# #                                       "color" : "black",
# #                                       "notes" : "GLODAP"
# #                                       },
# #                           "GLODAP_climatology" : {"var"  : var,
# #                                       "y0"    : 1980, 
# #                                       "y1"    : 2016,  
# #                                       "dir"   : f"{loc}" + f"/GLODAPv2.2016b.{var}.nc",
# #                                       "file"  : "GLODAP",
# #                                       "color" : "black",
# #                                       "notes" : "GLODAP"
# #                                       }
# #                                       }
# #           elif var in ['ntalk', 'talk','dissic', 'ndissic','o2']:
# #                       self.data_dict = { 
# #                           "GLODAP" : {"var"  : var,
# #                                       "y0"    : 1984, 
# #                                       "y1"    : 2021,  
# #                                       "dir"   : f"{loc}" + f"/GLODAP_v2_2023_{var}_1984-2021.csv",
# #                                       "file"  : "GLODAP",
# #                                       "color" : "black",
# #                                       "notes" : "GLODAP"
# #                                       },
# #                         "GLODAP_climatology" : {"var"  : var,
# #                                       "y0"    : 1984, 
# #                                       "y1"    : 2016,  
# #                                       "dir"   : f"{loc}" + f"/GLODAPv2.2016b.{var}.nc",
# #                                       "file"  : "GLODAP",
# #                                       "color" : "black",
# #                                       "notes" : "GLODAP"
# #                                       }}
# #           else:
# #                  self.data_dict = {}
                      
        
# #         elif all(['assimilation' in self.loc]):
# #               if var in ['uas','vas','tauu','tauv']:
# #                      realm = 'Amon'
# #               else:
# #                      realm = 'Omon'
# #               if 'CanOE' in self.loc:
# #                   color = "tab:blue"
# #               else:
# #                   color = "tab:green"

# #               self.data_dict = {"var"   : var,
# #                                       "y0"    : 1958, 
# #                                       "y1"    : 2023,  
# #                                       "dir"   : f"{loc}" + f"/{var}_{realm}_ensmebles_*_1x1_LE.nc",
# #                                       "file"  : "A1-assim-re",
# #                                       "color" : color,
# #                                       "notes" : "None"
# #                                       }       
                     
# #         elif all(['historical' in self.loc ]):

# #               if var in ['uas','vas','tauu','tauv']:
# #                      realm = 'Amon'
# #               else:
# #                      realm = 'Omon'
# #               if var == 'chlos':
# #                     y0 = 1982
# #               else:
# #                     y0 = 1950
# #               if 'CanOE' in self.loc:
# #                   color = "tab:purple"
# #               else:
# #                   color = "tab:red"
                     
# #               self.data_dict = { "var"   : var,
# #                                       "y0"    : y0, 
# #                                       "y1"    : 2025,  
# #                                       "dir"   : f"{loc}" + f"/{var}_{realm}_ensmebles_*_1x1_LE.nc",
# #                                       "file"  : "U1-simulation",
# #                                       "color" : color,
# #                                       "notes" : "None"
# #                                      }
              
# #         elif 'hindcast' in self.loc:


# #               if var == 'chlos':
# #                   self.data_dict = { "var"   : "chlos",
# #                                           "y0"    : 1991, 
# #                                           "y1"    : 2024,  
# #                                           "dir"   : f"{loc}",
# #                                           "file"  : "chlos_Omon",
# #                                           "color" : "tab:blue",
# #                                           "notes" : "None"
# #                                         }
# #               else:
# #                   self.data_dict = {}

# #         elif 'Control' in self.loc:       
# #               if var == 'chlos':
# #                   self.data_dict = {  "var"   : "chlos",
# #                                       "y0"    : 2801, 
# #                                       "y1"    : 3000,  
# #                                       "dir"   : f"{loc}"  + "/chlos_Omon_ensmebles_280101_300012_1x1_LE.nc", 
# #                                       "file"  : "Control_Omon",
# #                                       "color" : "grey",
# #                                       "notes" : "None"
# #                                      }
# #               else:
# #                   self.data_dict = {}
            


# #     def PrintLoc(self):
# #         print(self.loc)
        
                    



