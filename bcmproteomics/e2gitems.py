

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

psm_columns = ['AUC_UseFLAG', 'ActivationType', 'AddedBy', 'Charge', 'CreationTS',
               'DeltaCn', 'DeltaMassDa', 'DeltaMassPPM', 'DeltaScore', 'EXPRecNo',
               'EXPRunNo', 'EXPSearchNo', 'EXPTechRepNo', 'FirstScan', 'GeneCount',
               'GeneID', 'GeneList', 'HID', 'HIDCount', 'HIDList', 'IonInjectTime',
               'IonScore', 'IsolationInterference', 'LabelFLAG', 'MHDa', 'MSOrder',
               'MatchedIons', 'MissedCleavages', 'ModificationTS', 'Modifications',
               'PEP', 'PSMAmbiguity', 'PSM_IDG', 'PSM_UseFLAG', 'Peak_UseFLAG',
               'PeptRank', 'PrecursorArea', 'PrecursorArea_dstrAdj', 'ProteinCount',
               'ProteinGI', 'ProteinList', 'QuanInfo', 'QuanUsage', 'RTmin', 'Rank',
               'SearchEngineRank', 'Sequence', 'SequenceArea', 'SequenceModi',
               'SequenceModiCount', 'SpectrumFile', 'TaxonCount', 'TaxonID',
               'TaxonIDList', 'mzDa', 'oriFLAG', 'q_value'
]

tmt_columns = ['TMT_126', 'TMT_127_C', 'TMT_127_N', 'TMT_128_C', 'TMT_128_N', 'TMT_129_C',
            'TMT_129_N', 'TMT_130_C', 'TMT_130_N', 'TMT_131',
]
