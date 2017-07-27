% UGM
fprintf('Compiling UGM files...\n');
mex -IUGM/mex -outdir UGM/mex UGM/mex/UGM_makeEdgeVEC.c
mex -IUGM/mex -outdir UGM/mex UGM/mex/UGM_Decode_ICMC.c
mex -IUGM/mex -outdir UGM/mex UGM/mex/UGM_LogConfigurationPotentialC.c
mex -IUGM/mex -outdir UGM/mex UGM/mex/UGM_makeClampedPotentialsC.c
mex -IUGM/mex -outdir UGM/mex UGM/mex/UGM_Decode_GraphCutC.c
mex -IUGM/mex -outdir UGM/mex UGM/mex/UGM_Decode_AlphaExpansionC.c
mex -IUGM/mex -outdir UGM/mex UGM/mex/UGM_Decode_AlphaExpansionBetaShrinkC.c

% UGMep
fprintf('Compiling UGMep files...\n');
mex -IUGMep/mex -outdir UGMep/mex UGMep/mex/UGMep_Decode_ICMC.c
mex -IUGMep/mex -outdir UGMep/mex UGMep/mex/UGMep_EnergyC.c
mex -IUGMep/mex -outdir UGMep/mex UGMep/mex/UGMep_makeClampedEnergyC.c
mex -IUGMep/mex -outdir UGMep/mex UGMep/mex/UGMep_Decode_GraphCutC.c
mex -IUGMep/mex -outdir UGMep/mex UGMep/mex/UGMep_Decode_AlphaExpansionC.c
mex -IUGMep/mex -outdir UGMep/mex UGMep/mex/UGMep_Decode_ExpandShrinkC.c