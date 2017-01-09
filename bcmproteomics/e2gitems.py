

e2gcolumns = ['e2g_'+ c for c in ['GeneID', 'IDSet', 'IDGroup',
                                  'IDGroup_u2g', 'GPGroup', 'GPGroups_all', 'EXPLabelFLAG',
                                  'PSMs', 'PSMs_u2g', 'PeptidePrint',
                                  'PeptideCount', 'PeptideCount_u2g',
                                  'PeptideCount_S', 'PeptideCount_S_u2g',
                                  'nGPArea_Sum_cgpAdj', 'nGPArea_Sum_u2g',
                                  'nGPArea_Sum_u2g_all', 'nGPArea_Sum_max',
                                  'nGPArea_Sum_dstrAdj', #'GeneCapacity',  # need to add this back
                                  'n_iBAQ_dstrAdj',
]]

psm_columns = ['psm_AUC_UseFLAG', 'ActivationType', 'psm_AddedBy', 'Charge', 'psm_CreationTS',
               'DeltaCn', 'DeltaMassDa', 'DeltaMassPPM', 'DeltaScore', 'psm_EXPRecNo',
               'psm_EXPRunNo', 'psm_EXPSearchNo', 'FirstScan', 'psm_GeneCount',
               'psm_GeneID', 'psm_GeneList', 'psm_HID', 'psm_HIDCount', 'psm_HIDList',
               'IonInjectTime', 'IonScore', 'IsolationInterference', 'psm_LabelFLAG', 'MHDa',
               'MSOrder', 'MatchedIons', 'MissedCleavages', 'psm_ModificationTS', 'Modifications',
               'PEP', 'PSMAmbiguity', 'psm_PSM_IDG', 'psm_PSM_UseFLAG', 'psm_Peak_UseFLAG',
               'psm_PeptRank',
               'PrecursorArea', 'psm_PrecursorArea_split', 'psm_PrecursorArea_dstrAdj',
               'psm_ProteinCount',
               'psm_ProteinGI', 'psm_ProteinList', 'QuanInfo', 'QuanUsage', 'RTmin', 'Rank',
               'SearchEngineRank', 'Sequence', 'psm_SequenceArea', 'psm_SequenceModi',
               'psm_SequenceModiCount', 'SpectrumFile', 'psm_TaxonCount', 'psm_TaxonID',
               'psm_TaxonIDList', 'mzDa', 'psm_oriFLAG', 'q_value'
]

tmt_columns = ['TMT_126', 'TMT_127_C', 'TMT_127_N', 'TMT_128_C', 'TMT_128_N', 'TMT_129_C',
               'TMT_129_N', 'TMT_130_C', 'TMT_130_N', 'TMT_131',
]
