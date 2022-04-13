# -*- coding: utf8 -*-
"""
ClassifyMutPatt - Kazuki Nakamae 2020
Software pipeline for classifing mutation type predicted by InDelPhi
https://github.com/KazukiNakamae/ClassifyMutPatt
(c) 2020 Hiroshima University. All Rights Reserved.
"""

__version__ = "1.0.0"

from ClassifyMutPatt.ClassifyMutPattShared import *
# from MaChIAto.MaChIAtoClassifier import MaChIAtoClassifier, MaChIAtoPEClassifier
from ClassifyMutPatt.ClassifyMutPattFunctions import clean_files

def main():
    try:
        print("  \n===\ClassifyMutPatt/===")
        print("""
        -Software pipeline for classifing mutation type predicted by InDelPhi-
        """)
        print("Author: Kazuki Nakamae at Sakuma Tetsushi and Takashi Yamamoto lab, Hiroshima Univ, Japan")
        print("For support contact kazukinakamae@gmail.com")
        info("ClassifyMutPatt Version :" + __version__)

        parser = argparse.ArgumentParser(description="ClassifyMutPatt Parameters",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        # Arguments
        parser.add_argument("-i", "--input", type=str,  help="This parameter allows for the specification of the input file that contains the gene name and sequence.", required=True)
        parser.add_argument("-o", "--output", type=str,  help="This parameter allows for the specification of the output directory that contains the gene name and predicted mutation pattern (default: output).", default="output.csv", required=False)
        parser.add_argument("-s", "--sequence_column", type=str,  help="This parameter allows for the specification of the column name that contains sequences.", required=True)
        parser.add_argument("-g", "--gene_column", type=str,  help="This parameter allows for the specification of the column name that contains gene names. Each label should be unique.", required=True)
        parser.add_argument("-m", "--mismatch_column", type=str,  help="This parameter allows for the specification of the column names that contains mis-match information of target site.", required=True, action='append')
        parser.add_argument("-e", "--essential_gene_list", type=str,  help="This parameter allows for the specification of the file name contains essential gene names  (default: "").", default="", required=False)
        parser.add_argument("-c", "--cancer_gene_list", type=str,  help="This parameter allows for the specification of the file name contains oncogene names (default: "").", default="", required=False)
        args = parser.parse_args()

        # Read the input table
        info("Read input table...")
        input_pddf = pd.read_csv( args.input )
        input_pddf = input_pddf.rename(columns={args.sequence_column: "Evaluated_Sequence"})
        input_pddf = input_pddf.rename(columns={args.gene_column: "GeneName"})

        # Read the essential_gene_list
        info("Read essential_gene_list...")
        essential_pdsr = pd.read_csv( args.essential_gene_list ).transpose().iloc[0, :]

        # Read the cancer_gene_list
        info("Read cancer_gene_list...")
        oncogene_pdsr = pd.read_csv( args.cancer_gene_list, delimiter="\t" ).loc[ :, "GENE_NAME" ]

        # Make a hidden directory that is to contain the immediate files
        info("Prepare the system...")
        if not os.path.exists(".ClassifyMutPatt_temp"):
            if subprocess.check_call([ "mkdir", ".ClassifyMutPatt_temp" ]) != 0:
                raise ExecutionException( "Cannot make the .ClassifyMutPatt_temp directory." )
        
        # Threadpool Configuration
        os.environ['NUMEXPR_MAX_THREADS'] = '64'
        
        # Filter by mis-match
        info("Filter by the 1-2 bp mis-match number...")
        
        ### filtering
        # note: I could not use query() due to unclear error. This is an alternative method.
        # make mis-match table
        mm_pddf = input_pddf
        non_mm_pddf =  mm_pddf
        for col_name in args.mismatch_column:
            non_mm_pddf =  non_mm_pddf[  non_mm_pddf.loc[ :, col_name ] == 0 ]
        info( str(len(non_mm_pddf)) + " / " + str(len(input_pddf)) + " targets do not have a mis-match site." )
        info("Done")

        # Filter by essentiality
        info("Filter by the essentiality...")
        non_essential_pddf = non_mm_pddf
        for gene_index, gene_name in essential_pdsr.items():
            non_essential_pddf = non_essential_pddf.query( "~GeneName.str.contains('" + gene_name + "', na=False)", engine='python')
        info( str(len(non_essential_pddf)) + " / " + str(len(non_mm_pddf)) + " targets are not targeting an essential gene." )
        info("Done")

        # Filter by cancer-causing
        info("Filter by the cancer-causing...")
        non_oncogene_pddf = non_essential_pddf
        for gene_name in oncogene_pdsr.unique():
            non_oncogene_pddf = non_oncogene_pddf.query( "~GeneName.str.contains('" + gene_name + "', na=False)", engine='python')
        info( str(len(non_oncogene_pddf)) + " / " + str(len(non_essential_pddf)) + " targets are not oncogene." )
        info("Done")

        fltr_input_pddf = non_oncogene_pddf

        # Search targets that have the mutation pattrn predicted by InDelphi
        info("Search targets that have the outstanding mutation pattrn with InDelphi...")
        for index in np.arange( len( fltr_input_pddf ) ):

            temp_pdsr = fltr_input_pddf.iloc[ index, : ]
            pred_df, stats = hek293_inDelphi.predict( temp_pdsr[ "Evaluated_Sequence" ], 30 )
            stat_pddf = pd.DataFrame.from_dict( stats, orient="index" ).transpose()
            if index == 0:
                stat_hek293_pddf = stat_pddf
            else:
                stat_hek293_pddf = pd.concat( [stat_hek293_pddf, stat_pddf] )

        fltr_input_pddf = fltr_input_pddf.set_index( np.arange( len( fltr_input_pddf ) ) )
        stat_hek293_pddf = stat_hek293_pddf.set_index( np.arange( len( stat_hek293_pddf ) ) )
        complete_info_pddf = pd.concat( [fltr_input_pddf, stat_hek293_pddf], axis=1 )
        complete_info_pddf = complete_info_pddf.query( "`Precision` > 0.35" )
        info( str(len(complete_info_pddf)) + " / " + str(len(non_oncogene_pddf)) + " targets represents gRNAs where a small number of repair outcomes constitute nearly all outcomes (Precision > 0.35)." )
        
        # Classify targets by the mutation type
        info("Classify targets by the mutation type...")
        oneins_freq_3q = np.quantile( complete_info_pddf.loc[ :, "1-bp ins frequency"], 0.75 )
        # mmej_freq_3q = np.quantile( complete_info_pddf.loc[ :, "MH del frequency"], 0.75 )
        # nhej_freq_3q = np.quantile( complete_info_pddf.loc[ :, "MHless del frequency"], 0.75 )
        info( "The third quartile (Q3) of '1-bp ins frequency' is " + str( oneins_freq_3q ) )
        # info( "The third quartile (Q3) of 'MH del frequency' is " + str( mmej_freq_3q ) )
        # info( "The third quartile (Q3) of 'MHless del frequency' is " + str( nhej_freq_3q ) )
        # info( "Classify by the Q3s of all data set...")

        # Insertion
        info( "`1-bp ins frequency` > " + str( oneins_freq_3q ) )
        oneins_pddf = complete_info_pddf.query( "`1-bp ins frequency` > " + str( oneins_freq_3q ) )
        info( str(len(oneins_pddf)) + " / " + str(len(complete_info_pddf)) + " targets represent a donominant 1 bp insertion" )
        # A detailed classification of the 1bp ins by the base
        # TODO: 行インデックスでループ集計を実行し、メジャーな塩基の種類ごとに分ける。 
        ind_base_dict = { "A": 0, "T": 0, "C": 0, "G": 0 } # index for each base table
        for index in np.arange( len( oneins_pddf ) ):

            temp_pdsr = oneins_pddf.iloc[ index, : ]
            pred_df, stats = hek293_inDelphi.predict( temp_pdsr[ "Evaluated_Sequence" ], 30 )
            max_ins_ind = pred_df.query( "Category.str.contains('ins', na=False)", engine='python' ).loc[ :, "Predicted frequency" ].idxmax()
            pred_max_pdsr = pred_df.iloc[ max_ins_ind, : ]
            ind_base_dict[ pred_max_pdsr[ "Inserted Bases" ] ] = ind_base_dict[ pred_max_pdsr[ "Inserted Bases" ] ] + 1
            temp_pddf = pd.DataFrame( pd.concat( [temp_pdsr, pred_max_pdsr], axis=0 ) ).transpose()

            if ( temp_pddf[ "Inserted Bases" ] == "G" ).bool():

                if ind_base_dict[ temp_pddf.loc[ 0, "Inserted Bases" ] ] == 1:

                    insG_hek293_pddf = temp_pddf

                else:

                    insG_hek293_pddf = pd.concat( [ insG_hek293_pddf, temp_pddf ], axis=0 )

            elif ( temp_pddf[ "Inserted Bases" ] == "C" ).bool():

                if ind_base_dict[ temp_pddf.loc[ 0, "Inserted Bases" ] ] == 1:

                    insC_hek293_pddf = temp_pddf

                else:
                    
                    insC_hek293_pddf = pd.concat( [ insC_hek293_pddf, temp_pddf ], axis=0 )

            elif ( temp_pddf[ "Inserted Bases" ] == "A" ).bool():

                if ind_base_dict[ temp_pddf.loc[ 0, "Inserted Bases" ] ] == 1:

                    insA_hek293_pddf = temp_pddf

                else:
                    
                    insA_hek293_pddf = pd.concat( [ insA_hek293_pddf, temp_pddf ], axis=0 )

            elif ( temp_pddf[ "Inserted Bases" ] == "T" ).bool():

                if ind_base_dict[ temp_pddf.loc[ 0, "Inserted Bases" ] ] == 1:

                    insT_hek293_pddf = temp_pddf

                else:
                    
                    insT_hek293_pddf = pd.concat( [ insT_hek293_pddf, temp_pddf ], axis=0 )

        if 'insA_hek293_pddf' in locals():
            info( str(len(insA_hek293_pddf)) + " / " + str(len(oneins_pddf)) + " targets represent a donominat 1 bp A insertion" )
        if 'insT_hek293_pddf' in locals():
            info( str(len(insT_hek293_pddf)) + " / " + str(len(oneins_pddf)) + " targets represent a donominat 1 bp T insertion" )
        if 'insC_hek293_pddf' in locals():
            info( str(len(insC_hek293_pddf)) + " / " + str(len(oneins_pddf)) + " targets represent a donominat 1 bp C insertion" )
        if 'insG_hek293_pddf' in locals():
            info( str(len(insG_hek293_pddf)) + " / " + str(len(oneins_pddf)) + " targets represent a donominat 1 bp G insertion" )
        
        # MMEJ
        info( "`1-bp ins frequency` <= " + str( oneins_freq_3q ) + " and  `MH del frequency` > `MHless del frequency`" )
        mmej_pddf = complete_info_pddf.query( "`1-bp ins frequency` <= " + str( oneins_freq_3q ) + " and  `MH del frequency` > `MHless del frequency`" )
        info( str(len(mmej_pddf)) + " / " + str(len(complete_info_pddf)) + " targets" )
        
        # MMEJ frame-in
        info( "`Frame +0 frequency` > `Frame +1 frequency` and `Frame +0 frequency` > `Frame +2 frequency`" )
        mmej_framein_pddf = mmej_pddf.query( "`Frame +0 frequency` > `Frame +1 frequency` and `Frame +0 frequency` > `Frame +2 frequency`" )
        info( str(len(mmej_framein_pddf)) + " / " + str(len(mmej_pddf)) + " targets represent a dominant MMEJ-assisted frame-in" )
        # MMEJ frameshift +1
        info( "`Frame +1 frequency` > `Frame +0 frequency` and `Frame +1 frequency` > `Frame +2 frequency`" )
        mmej_framep1_pddf = mmej_pddf.query( "`Frame +1 frequency` > `Frame +0 frequency` and `Frame +1 frequency` > `Frame +2 frequency`" )
        info( str(len(mmej_framep1_pddf)) + " / " + str(len(mmej_pddf)) + " targets represent a dominant MMEJ-assisted frameshift (+1)" )
        # MMEJ frameshift +2
        info( "`Frame +2 frequency` > `Frame +0 frequency` and `Frame +2 frequency` > `Frame +1 frequency`" )
        mmej_framep2_pddf = mmej_pddf.query( "`Frame +2 frequency` > `Frame +0 frequency` and `Frame +2 frequency` > `Frame +1 frequency`" )
        info( str(len(mmej_framep2_pddf)) + " / " + str(len(mmej_pddf)) + " targets represent a dominant MMEJ-assisted frameshift (+2)" )

        # NHEJ
        info( "`1-bp ins frequency` <= " + str( oneins_freq_3q ) + " and  `MHless del frequency` > `MH del frequency`" )
        nhej_pddf = complete_info_pddf.query( "`1-bp ins frequency` <= " + str( oneins_freq_3q ) + " and  `MHless del frequency` > `MH del frequency`" )
        info( str(len(nhej_pddf)) + " / " + str(len(complete_info_pddf)) + " targets" )
        
        # nhej frame-in
        info( "`Frame +0 frequency` > `Frame +1 frequency` and `Frame +0 frequency` > `Frame +2 frequency`" )
        nhej_framein_pddf = nhej_pddf.query( "`Frame +0 frequency` > `Frame +1 frequency` and `Frame +0 frequency` > `Frame +2 frequency`" )
        info( str(len(nhej_framein_pddf)) + " / " + str(len(nhej_pddf)) + " targets represent a dominant NHEJ-assisted frame-in" )
        # nhej frameshift +1
        info( "`Frame +1 frequency` > `Frame +0 frequency` and `Frame +1 frequency` > `Frame +2 frequency`" )
        nhej_framep1_pddf = nhej_pddf.query( "`Frame +1 frequency` > `Frame +0 frequency` and `Frame +1 frequency` > `Frame +2 frequency`" )
        info( str(len(nhej_framep1_pddf)) + " / " + str(len(nhej_pddf)) + " targets represent a dominant NHEJ-assisted frameshift (+1)" )
        # nhej frameshift +2
        info( "`Frame +2 frequency` > `Frame +0 frequency` and `Frame +2 frequency` > `Frame +1 frequency`" )
        nhej_framep2_pddf = nhej_pddf.query( "`Frame +2 frequency` > `Frame +0 frequency` and `Frame +2 frequency` > `Frame +1 frequency`" )
        info( str(len(nhej_framep2_pddf)) + " / " + str(len(nhej_pddf)) + " targets represent a dominant NHEJ-assisted frameshift (+2)" )

        info("Sort...")
        # 1 bp ins A
        sorted_insA_hek293_pddf = insA_hek293_pddf.sort_values("Predicted frequency", ascending=False)

        # 1 bp ins T
        sorted_insT_hek293_pddf = insT_hek293_pddf.sort_values("Predicted frequency", ascending=False)

        # 1 bp ins C
        sorted_insC_hek293_pddf = insC_hek293_pddf.sort_values("Predicted frequency", ascending=False)

        # MMEJ / frame-in
        mmej_framein_pddf.insert( len(mmej_framein_pddf.columns) , "MHdel/MHlessdel", mmej_framein_pddf.loc[ :, "MH del frequency" ] / mmej_framein_pddf.loc[ :, "MHless del frequency" ] )
        sorted_mmej_framein_hek293_pddf = mmej_framein_pddf.sort_values("MHdel/MHlessdel", ascending=False)
        
        # MMEJ / frameshift +1
        mmej_framep1_pddf.insert( len(mmej_framep1_pddf.columns) , "MHdel/MHlessdel", mmej_framep1_pddf.loc[ :, "MH del frequency" ] / mmej_framep1_pddf.loc[ :, "MHless del frequency" ] )
        sorted_mmej_framep1_hek293_pddf = mmej_framep1_pddf.sort_values("MHdel/MHlessdel", ascending=False)
        
        # MMEJ / frameshift +2
        mmej_framep2_pddf.insert( len(mmej_framep2_pddf.columns) , "MHdel/MHlessdel", mmej_framep2_pddf.loc[ :, "MH del frequency" ] / mmej_framep2_pddf.loc[ :, "MHless del frequency" ] )
        sorted_mmej_framep2_hek293_pddf = mmej_framep2_pddf.sort_values("MHdel/MHlessdel", ascending=False)
        
        # NHEJ / frame-in
        nhej_framein_pddf.insert( len(nhej_framein_pddf.columns) , "MHlessdel/MHdel", nhej_framein_pddf.loc[ :, "MHless del frequency" ] / nhej_framein_pddf.loc[ :, "MH del frequency" ] )
        sorted_nhej_framein_hek293_pddf = nhej_framein_pddf.sort_values("MHlessdel/MHdel", ascending=False)
        
        # NHEJ / frameshift +1
        nhej_framep1_pddf.insert( len(nhej_framep1_pddf.columns) , "MHlessdel/MHdel", nhej_framep1_pddf.loc[ :, "MHless del frequency" ] / nhej_framep1_pddf.loc[ :, "MH del frequency" ] )
        sorted_nhej_framep1_hek293_pddf = nhej_framep1_pddf.sort_values("MHlessdel/MHdel", ascending=False)
        
        # NHEJ / frameshift +2
        nhej_framep2_pddf.insert( len(nhej_framep2_pddf.columns) , "MHlessdel/MHdel", nhej_framep2_pddf.loc[ :, "MHless del frequency" ] / nhej_framep2_pddf.loc[ :, "MH del frequency" ] )
        sorted_nhej_framep2_hek293_pddf = nhej_framep2_pddf.sort_values("MHlessdel/MHdel", ascending=False)

        info("Evaluating...")
        def ProfileMutAcrossCellline( row ):

            def LabelMut( pred_row ):
                if ( pred_row["Category"] is "del" ) and ( pred_row["Genotype position"] is not "e" ):
                    return( "LEN" + str( pred_row["Length"] ) + "R" + str( pred_row["Genotype position"] ) + pred_row["Category"] )
                elif ( pred_row["Category"] is "del" ) and ( pred_row["Genotype position"] is "e" ):
                    return( "LEN" + str( pred_row["Length"] ) + "NHEJ" + pred_row["Category"] )
                elif ( pred_row["Category"] is "ins" ):
                    return( str( pred_row["Inserted Bases"] ) + pred_row["Category"] )
                else:
                    raise UnexpectedException('Unknown pattrn is detected.')
            
            # HEK293
            hek293_pred_df, hek293_stats = hek293_inDelphi.predict( row[ "Evaluated_Sequence" ], 30 )
            named_hek293_pred_df = hek293_pred_df.set_index( hek293_pred_df.apply( LabelMut, axis=1 ) )
            named_hek293_pred_df = named_hek293_pred_df.rename( columns={ "Predicted frequency": "HEK293 Predicted frequency" } )
            # U2OS
            u2os_pred_df, u2os_stats = u2os_inDelphi.predict( row[ "Evaluated_Sequence" ], 30 )
            named_u2os_pred_df = u2os_pred_df.set_index( u2os_pred_df.apply( LabelMut, axis=1 ) )
            named_u2os_pred_df = named_u2os_pred_df.rename( columns={ "Predicted frequency": "U2OS Predicted frequency" } )
            # HCT116
            hct116_pred_df, hct116_stats = hct116_inDelphi.predict( row[ "Evaluated_Sequence" ], 30 )
            named_hct116_pred_df = hct116_pred_df.set_index( hct116_pred_df.apply( LabelMut, axis=1 ) )
            named_hct116_pred_df = named_hct116_pred_df.rename( columns={ "Predicted frequency": "HCT116 Predicted frequency" } )
            # K562
            k562_pred_df, k562_stats = k562_inDelphi.predict( row[ "Evaluated_Sequence" ], 30 )
            named_k562_pred_df = k562_pred_df.set_index( k562_pred_df.apply( LabelMut, axis=1 ) )
            named_k562_pred_df = named_k562_pred_df.rename( columns={ "Predicted frequency": "K562 Predicted frequency" } )
            
            # Calculate (KL Divergence P||Q) = SUM(P*log(P/Q))
            ### Test
            # P = torch.Tensor([0.36, 0.48, 0.16, 0])
            # Q = torch.Tensor([0.333, 0.333, 0.333, 0.1])
            # kldiv(Q.log(), P)
            ### Test END
            # Merge data according to the mutation predicted by hek293 model
            merged__hek293_pred_df = named_hek293_pred_df.loc[ :, [ "HEK293 Predicted frequency" ] ].join( named_u2os_pred_df.loc[ :, [ "U2OS Predicted frequency" ] ] )
            merged__hek293_pred_df = merged__hek293_pred_df.join( named_hct116_pred_df.loc[ :, [ "HCT116 Predicted frequency" ] ] )
            merged__hek293_pred_df = merged__hek293_pred_df.join( named_k562_pred_df.loc[ :, [ "K562 Predicted frequency" ] ] )
            torch.Tensor( merged__hek293_pred_df.loc[ :, [ "HEK293 Predicted frequency" ] ].transpose().values )
            # Calcutate D_kl(HEK293 mutation distribution || Other cell mutation distribution)
            hek293mut_hek293freq = torch.Tensor( merged__hek293_pred_df.loc[ :, [ "HEK293 Predicted frequency" ] ].values )
            hek293mut_u2osfreq = torch.Tensor( merged__hek293_pred_df.loc[ :, [ "U2OS Predicted frequency" ] ].values )
            hek293mut_hct116freq = torch.Tensor( merged__hek293_pred_df.loc[ :, [ "HCT116 Predicted frequency" ] ].values )
            hek293mut_k562freq = torch.Tensor( merged__hek293_pred_df.loc[ :, [ "K562 Predicted frequency" ] ].values )
            D_kl_hek293_hek293 = kldiv(hek293mut_hek293freq.log(), hek293mut_hek293freq) # Dkl(HEK293 || HEK293)
            D_kl_hek293_u2os = kldiv(hek293mut_u2osfreq.log(), hek293mut_hek293freq) # Dkl(HEK293 || U2OS)
            D_kl_hek293_hct116 = kldiv(hek293mut_hct116freq.log(), hek293mut_hek293freq) # Dkl(HEK293 || HCT116)
            D_kl_hek293_k562 = kldiv(hek293mut_k562freq.log(), hek293mut_hek293freq) # Dkl(HEK293 || K562)
            kl_list = [ D_kl_hek293_hek293.item(), D_kl_hek293_u2os.item(), D_kl_hek293_hct116.item(), D_kl_hek293_k562.item() ]
            # Calcutate mean of KL divergence. The applied calculatation is arithmetic mean.
            # Reference: https://doi.org/10.1177%2F1550147717747848
            kl_ave = sum( kl_list ) / len( kl_list )
            # Return values
            kl_list.append( kl_ave )
            return( pd.Series( kl_list, index =['Dkl_hek293_hek293', 'Dkl_hek293_u2os', 'Dkl_hek293_hct116', 'Dkl_hek293_k562', 'Dkl_ave'] ) )

        # import pdb; pdb.set_trace() #---BREAKPOINT----------------------------------------------------------------------------------------------#
        # Calculate KL divergence
        evaluated_insA_hek293_pddf = pd.concat( [ sorted_insA_hek293_pddf.head( 10 ), sorted_insA_hek293_pddf.head( 10 ).apply( ProfileMutAcrossCellline, axis=1 ) ], axis=1 )
        evaluated_insT_hek293_pddf = pd.concat( [ sorted_insT_hek293_pddf.head( 10 ), sorted_insT_hek293_pddf.head( 10 ).apply( ProfileMutAcrossCellline, axis=1 ) ], axis=1 )
        evaluated_insC_hek293_pddf = pd.concat( [ sorted_insC_hek293_pddf.head( 10 ), sorted_insC_hek293_pddf.head( 10 ).apply( ProfileMutAcrossCellline, axis=1 ) ], axis=1 )
        evaluated_mmej_framein_hek293_pddf = pd.concat( [ sorted_mmej_framein_hek293_pddf.head( 10 ), sorted_mmej_framein_hek293_pddf.head( 10 ).apply( ProfileMutAcrossCellline, axis=1 ) ], axis=1 )
        evaluated_mmej_framep1_hek293_pddf = pd.concat( [ sorted_mmej_framep1_hek293_pddf.head( 10 ), sorted_mmej_framep1_hek293_pddf.head( 10 ).apply( ProfileMutAcrossCellline, axis=1 ) ], axis=1 )
        evaluated_mmej_framep2_hek293_pddf = pd.concat( [ sorted_mmej_framep2_hek293_pddf.head( 10 ), sorted_mmej_framep2_hek293_pddf.head( 10 ).apply( ProfileMutAcrossCellline, axis=1 ) ], axis=1 )
        # evaluated_nhej_framein_hek293_pddf = pd.concat( [ sorted_nhej_framein_hek293_pddf.head( 10 ), sorted_nhej_framein_hek293_pddf.head( 10 ).apply( ProfileMutAcrossCellline, axis=1 ) ], axis=1 )
        evaluated_nhej_framep1_hek293_pddf = pd.concat( [ sorted_nhej_framep1_hek293_pddf.head( 5 ), sorted_nhej_framep1_hek293_pddf.head( 5 ).apply( ProfileMutAcrossCellline, axis=1 ) ], axis=1 )
        # evaluated_nhej_framep2_hek293_pddf = pd.concat( [ sorted_nhej_framep2_hek293_pddf.head( 10 ), sorted_nhej_framep2_hek293_pddf.head( 10 ).apply( ProfileMutAcrossCellline, axis=1 ) ], axis=1 )
        # Summarize result
        comp_target_pddf = pd.concat( [ \
            evaluated_insA_hek293_pddf[ evaluated_insA_hek293_pddf.Dkl_ave == evaluated_insA_hek293_pddf.Dkl_ave.max() ].assign( Value_name="insA_Dkl_max" )\
            , evaluated_insA_hek293_pddf[ evaluated_insA_hek293_pddf.Dkl_ave == evaluated_insA_hek293_pddf.Dkl_ave.min() ].assign( Value_name="insA_Dkl_min" )\
            , evaluated_insT_hek293_pddf[ evaluated_insT_hek293_pddf.Dkl_ave == evaluated_insT_hek293_pddf.Dkl_ave.max() ].assign( Value_name="insT_Dkl_max" )\
            , evaluated_insT_hek293_pddf[ evaluated_insT_hek293_pddf.Dkl_ave == evaluated_insT_hek293_pddf.Dkl_ave.min() ].assign( Value_name="insT_Dkl_min" )\
            , evaluated_insC_hek293_pddf[ evaluated_insC_hek293_pddf.Dkl_ave == evaluated_insC_hek293_pddf.Dkl_ave.max() ].assign( Value_name="insC_Dkl_max" )\
            , evaluated_insC_hek293_pddf[ evaluated_insC_hek293_pddf.Dkl_ave == evaluated_insC_hek293_pddf.Dkl_ave.min() ].assign( Value_name="insC_Dkl_min" )\
            , evaluated_mmej_framein_hek293_pddf[ evaluated_mmej_framein_hek293_pddf.Dkl_ave == evaluated_mmej_framein_hek293_pddf.Dkl_ave.max() ].assign( Value_name="mmej_framein_Dkl_max" )\
            , evaluated_mmej_framein_hek293_pddf[ evaluated_mmej_framein_hek293_pddf.Dkl_ave == evaluated_mmej_framein_hek293_pddf.Dkl_ave.min() ].assign( Value_name="mmej_framein_Dkl_min" )\
            , evaluated_mmej_framep1_hek293_pddf[ evaluated_mmej_framep1_hek293_pddf.Dkl_ave == evaluated_mmej_framep1_hek293_pddf.Dkl_ave.max() ].assign( Value_name="mmej_framep1_Dkl_max" )\
            , evaluated_mmej_framep1_hek293_pddf[ evaluated_mmej_framep1_hek293_pddf.Dkl_ave == evaluated_mmej_framep1_hek293_pddf.Dkl_ave.min() ].assign( Value_name="mmej_framep1_Dkl_min" )\
            , evaluated_mmej_framep2_hek293_pddf[ evaluated_mmej_framep2_hek293_pddf.Dkl_ave == evaluated_mmej_framep2_hek293_pddf.Dkl_ave.max() ].assign( Value_name="mmej_framep2_Dkl_max" )\
            , evaluated_mmej_framep2_hek293_pddf[ evaluated_mmej_framep2_hek293_pddf.Dkl_ave == evaluated_mmej_framep2_hek293_pddf.Dkl_ave.min() ].assign( Value_name="mmej_framep2_Dkl_min" )\
            # , evaluated_nhej_framein_hek293_pddf[ evaluated_nhej_framein_hek293_pddf.Dkl_ave == evaluated_nhej_framein_hek293_pddf.Dkl_ave.max() ].assign( Value_name="nhej_framein_Dkl_max" )\
            # , evaluated_nhej_framein_hek293_pddf[ evaluated_nhej_framein_hek293_pddf.Dkl_ave == evaluated_nhej_framein_hek293_pddf.Dkl_ave.min() ].assign( Value_name="nhej_framein_Dkl_min" )\
            , evaluated_nhej_framep1_hek293_pddf[ evaluated_nhej_framep1_hek293_pddf.Dkl_ave == evaluated_nhej_framep1_hek293_pddf.Dkl_ave.max() ].assign( Value_name="nhej_framep1_Dkl_max" )\
            , evaluated_nhej_framep1_hek293_pddf[ evaluated_nhej_framep1_hek293_pddf.Dkl_ave == evaluated_nhej_framep1_hek293_pddf.Dkl_ave.min() ].assign( Value_name="nhej_framep1_Dkl_min" )\
            # , evaluated_nhej_framep2_hek293_pddf[ evaluated_nhej_framep2_hek293_pddf.Dkl_ave == evaluated_nhej_framep2_hek293_pddf.Dkl_ave.max() ].assign( Value_name="nhej_framep2_Dkl_max" )\
            # , evaluated_nhej_framep2_hek293_pddf[ evaluated_nhej_framep2_hek293_pddf.Dkl_ave == evaluated_nhej_framep2_hek293_pddf.Dkl_ave.min() ].assign( Value_name="nhej_framep2_Dkl_min" )\
            ], axis=0 )
        
        # Save result
        info("Saving...")
        if subprocess.check_call([ "mkdir", args.output ]) != 0:
            raise ExecutionException( "Cannot make the output directory." )
        comp_target_pddf.to_csv( path_or_buf=args.output + "/selected_target_list.csv", sep=",", index=False )
        evaluated_insA_hek293_pddf.to_csv( path_or_buf=args.output + "/insA_target_list.csv", sep=",", index=False )
        evaluated_insT_hek293_pddf.to_csv( path_or_buf=args.output + "/insT_target_list.csv", sep=",", index=False )
        evaluated_insC_hek293_pddf.to_csv( path_or_buf=args.output + "/insC_target_list.csv", sep=",", index=False )
        evaluated_mmej_framein_hek293_pddf.to_csv( path_or_buf=args.output + "/mmej_framein_target_list.csv", sep=",", index=False )
        evaluated_mmej_framep1_hek293_pddf.to_csv( path_or_buf=args.output + "/mmej_framep1_target_list.csv", sep=",", index=False )
        evaluated_mmej_framep2_hek293_pddf.to_csv( path_or_buf=args.output + "/mmej_framep2_target_list.csv", sep=",", index=False )
        # evaluated_nhej_framein_hek293_pddf.to_csv( path_or_buf=args.output + "/nhej_framein_target_list.csv", sep=",", index=False )
        evaluated_nhej_framep1_hek293_pddf.to_csv( path_or_buf=args.output + "/nhej_framep1_target_list.csv", sep=",", index=False )
        # evaluated_nhej_framep2_hek293_pddf.to_csv( path_or_buf=args.output + "/nhej_framep2_target_list.csv", sep=",", index=False )

        # Remove the immediate files
        info("The process is ending...")
        if clean_files():
            info("Done.")

        sys.exit(0)

    except Exception as e:
        # import pdb; pdb.set_trace() #---BREAKPOINT----------------------------------------------------------------------------------------------#
        error("Error. Please check your input.\n\nERROR: %s" % e)
        # clean_files()
        sys.exit(-1)
    except InputException as e:
        #traceback.print_exc(file = sys.stdout)
        error("Input error. Please check your input.\n\nERROR: %s" % e)
        print("""
        usage:
        python ExtendSeq.py -i <INPUT> -o <OUTPUT> -n <NAME COLUMN> -s <SEQUENCE COLUMN> -r <REFERENCE COLUMN> -l <EXTEND LEFT LENGTH> -r <EXTENDED RIGHT LENGTH>
        """)
        clean_files()
        sys.exit(1)
    except ExecutionException as e:
        #traceback.print_exc(file = sys.stdout)
        error("Execution error. Please check your input.\n\nERROR: %s" % e)
        clean_files()
        sys.exit(1)
    except UnexpectedException as e:
        #traceback.print_exc(file = sys.stdout)
        error("Unexpected error. Please get touch in the author.[Kazuki Nakamae : kazukinakamae@gmail.com]\n\nERROR: %s" % e)
        # clean_files()
        sys.exit(1)
