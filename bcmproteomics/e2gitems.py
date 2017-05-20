

# e2gcolumns = ['e2g_'+ c for c in ['GeneID', 'IDSet', 'IDGroup', 'TaxonID',
#                                   'IDGroup_u2g', 'GPGroup', 'GPGroups_all', 'EXPLabelFLAG',
#                                   'PSMs', 'PSMs_u2g', 'PeptidePrint',
#                                   'PeptideCount', 'PeptideCount_u2g',
#                                   'PeptideCount_S', 'PeptideCount_S_u2g',
#                                   'nGPArea_Sum_cgpAdj', 'nGPArea_Sum_u2g',
#                                   'nGPArea_Sum_u2g_all', 'nGPArea_Sum_max',
#                                   'nGPArea_Sum_dstrAdj', 'GeneCapacity',
#                                   'n_iBAQ_dstrAdj',
# ]]


e2gcolumns = ['e2g_' + c for c in (
    'GeneID',
    'TaxonID',
    'IDSet',
    'IDGroup',
    'IDGroup_u2g',
    'GPGroup',
    'GPGroups_all',
    'EXPLabelFLAG',
    'PSMs',
    'PSMs_u2g',
    'PeptidePrint',
    'PeptideCount',
    'PeptideCount_u2g',
    'PeptideCount_S',
    'PeptideCount_S_u2g',
    'AreaSum_gpcAdj',
    'AreaSum_u2g_0',
    'AreaSum_u2g_all',
    'AreaSum_max',
    'AreaSum_dstrAdj',
    'AreaSum_razor',
    'GeneCapacity',
    'iBAQ_dstrAdj',
)]


psm_columns = ['psm_' + c for c in (
    'EXPRecNo',
    'EXPRunNo',
    'EXPSearchNo',
    'Sequence',
    'PSMAmbiguity',
    'Modifications',
    'ActivationType',
    'DeltaScore',
    'DeltaCn',
    'Rank',
    'SearchEngineRank',
    'PrecursorArea',
    'q_value',
    'PEP',
    'IonScore',
    'MissedCleavages',
    'IonInjectTime',
    'Charge',
    'mzDa',
    'MHDa',
    'DeltaMassDa',
    'DeltaMassPPM',
    'RTmin',
    'FirstScan',
    'MSOrder',
    'QuanInfo',
    'QuanUsage',
    'SpectrumFile',
    'AddedBy',
    'CreationTS',
    'ModificationTS',
    'oriFLAG',
    'LabelFLAG',
    'PSM_IDG',
    'SequenceModi',
    'SequenceModiCount',
    'PSM_UseFLAG',
    'AUC_UseFLAG',
    'Peak_UseFLAG',
    'PrecursorArea_split',
    'SequenceArea',
    'PrecursorArea_dstrAdj',
    'RazorArea',
    'PeptRank',
    'GeneID',
    'GeneIDs_All',
    'GeneIDCount_All',
    'ProteinGIs',
    'ProteinGIs_All',
    'ProteinGICount_All',
    'ProteinRefs',
    'ProteinRefs_All',
    'ProteinRefCount_All',
    'HIDs',
    'HIDCount_All',
    'TaxonID',
    'TaxonIDs_All',
    'TaxonIDCount_All'
)]


tmt_columns = ['psm_' + x for x in ['TMT_126', 'TMT_127_C', 'TMT_127_N', 'TMT_128_C', 'TMT_128_N', 'TMT_129_C',
               'TMT_129_N', 'TMT_130_C', 'TMT_130_N', 'TMT_131',]
]
