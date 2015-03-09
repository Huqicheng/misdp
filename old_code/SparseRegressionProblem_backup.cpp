//
//  SparseRegressionProblem.cpp
//  SparseRegression
//
//  Created by Qi on 5/27/14.
//  Copyright (c) 2014 Qi Zhang. All rights reserved.
//

/*
#define debugMode_display_stage 1
#define debugMode_positionTitle 1

#define debugMode_display_cut_info_summary 1
#define pause_20140604_0744 0
#define pauseUpperGrad 1
#define debugMode_display_all_cuts 1

#define debugMode 1
#define debugMode_usercut 1
#define debugMode_usercut_gradient 1
#define debugMode_display_callback_info 1
*/

//#define debugMode_display_callback_sol 1
//#define debugMode_display_sol 1
//#define debugMode_display_round1 1
//#define debugMode_display_cut 1
//#define debugMode_display_pause 1


//#define pause_total_grad 1
//#define pause_total_grad2 1
//#define pause_total_grad3 1
//#define pause_lower_facet 1
//#define pause_exp_grad1 1
//
//#define pause_initial_cons 1
//#define pause_total_grad4 1
//
//
//



//#define title_everywhere 1
//#define debug_20140911_1pm_findSegmentationFault 1






// add on 10/10/2014
#define test_nsd_cut 0
#define test_nsd_cut_level2 0
#define pause_at_incumbent 0
//#define test_my_node 0

#define test_upgrad 0

#define test_psd_cut 0


#define test_adap_psd_cut 1
#define test_adap_psd_cut_level2 1
#define test_adap_psd_cut_level3 1
#define test_adap_psd_cut_level4 1
#define test_adap_psd_cut_level5 1
#define test_adap_nsd_cut 0
#define test_adap_nsd_cut_l2 0


#define debugMode_display_cut_info_summary 1

//const int test_array [ 8 ] = {1,1,0,1,1,0,0,0};

#include "SparseRegressionProblem.h"

using namespace std;

void SparseRegressionProblem :: generateCuts_lowercuts( IloEnv &env, IloConstraintArray &cuts2BeAdded, DataModel &data, TuneParameters &tune, IloBoolVarArray &z, IloNumVar &g, IloNumArray &valz, IloNum valg)
{
    
#if title_everywhere
    cout << " SparseRegressionProblem :: generateCuts_lowercuts tag (f4w98f) " << endl;
#endif
    int n = data.getCols();
//    double cumulative_violated_amount = 0;
    IloRange cut;
    IloExpr  rhs(env);
    IloNum   valrhs = 0;
#if pause_20140604_0744
    IloNum test_rhs = 0;
#endif

    bitset<_CONST_INT_MAX_FACTOR_SIZE> chainsetFront;
    bitset<_CONST_INT_MAX_FACTOR_SIZE> chainsetBack;
    
    vector<int> sorted_index( Qi_Algorithm::sort_IloNumArray_sub_indexes_decreasing( valz, 0, n ) );
    
#if debugMode_display_sort
    {
        for (int temp = 0; temp < sorted_index.size(); temp++)
            cout << sorted_index.at(temp) << "\t" ;
        cout << endl;
    }
#endif
    
    try {
        
        rhs += data.logDetSubmatrix_with_y(chainsetFront);
        valrhs += data.logDetSubmatrix_with_y(chainsetFront);
#if pause_20140604_0744
        test_rhs += data.logDetSubmatrix_with_y( chainsetFront );
#endif

        for ( int j : sorted_index ) {
            chainsetFront.set( j );
            rhs += ( data.logDetSubmatrix_with_y(chainsetFront) - data.logDetSubmatrix_with_y(chainsetBack) ) * z[j];
            valrhs += ( data.logDetSubmatrix_with_y(chainsetFront) - data.logDetSubmatrix_with_y(chainsetBack) ) * valz[j];
#if pause_20140604_0744
            test_rhs += ( data.logDetSubmatrix_with_y(chainsetFront) - data.logDetSubmatrix_with_y(chainsetBack) ) * test_array[j];;
#endif
            chainsetBack = chainsetFront;
        };
        cut = ( g - rhs >= 0 );
        tune.cumulative_violated_amount +=  abs(valg - valrhs) ;
        
        
#if pause_20140604_0744
        cout <<" - - {tag iojqf0q234}- - - { < 19.06319} test_rhs = " << test_rhs << endl;
        if ( test_rhs > 19.06319 + 0.01 ){
            cout << " - something wrong here ( tag f35v57nh4w5byw45 )" << endl ;
            int temp;
            cin >> temp;
        }
#endif

        
#if debugMode_display_cutoff
        cout <<  ", Lcut cut off, g - rhs = " << valg - valrhs << endl;
#endif
        
#if debugMode_display_cut || debugMode_display_all_cuts
        cout << "\n %% lower facet cut (f34v56b) : " << cut << endl;
//        cout << " @@ currently, valg = " << valg << ", rhs = " << valrhs << endl << endl;
#endif
        
#if debugMode_usercut_gradient
        {
            cout <<"tag 872439r81yr92874" << endl;
            cout << " valz = " << valz << endl;
            cout << " valrhs = " << valrhs << ", valg = " << valg << endl;
            int tempp;
            cin >> tempp;
        }
#endif

        if (valg < valrhs - tune.DOUBLE_EPSILON_ADD_USER_CUT_THRESHOLD - 0.001){
#if pause_lower_facet
            cout << "lower facet cut = " << endl;
            cout << cut << endl;
#endif
            cuts2BeAdded.add(cut);
            tune.num_sbm_lower_cut_on_g++;
        }
        
#if debugMode_display_callback_sol
        cout << " valz = " << valz;
        cout << ", valg = " << valg << endl;
#endif

        
#if debugMode_display_pause
        {
            cout << " - have a look on the cut above ( 342cu0352 )" << endl ;
            int temp;
            cin >> temp;
        }
#endif
        
        tune.any_cut_added=true;
        tune.num_cut_added_this_round++;
        
        rhs.end();
        
    } catch (...) {
        rhs.end();
        
        throw;
    }

};



//               ###    ########     ###    ########  ######## #### ##     ## ########    ########  #######  ########    ###    ##
//              ## ##   ##     ##   ## ##   ##     ##    ##     ##  ##     ## ##             ##    ##     ##    ##      ## ##   ##
//             ##   ##  ##     ##  ##   ##  ##     ##    ##     ##  ##     ## ##             ##    ##     ##    ##     ##   ##  ##
//            ##     ## ##     ## ##     ## ########     ##     ##  ##     ## ######         ##    ##     ##    ##    ##     ## ##
//            ######### ##     ## ######### ##           ##     ##   ##   ##  ##             ##    ##     ##    ##    ######### ##
//            ##     ## ##     ## ##     ## ##           ##     ##    ## ##   ##             ##    ##     ##    ##    ##     ## ##
//            ##     ## ########  ##     ## ##           ##    ####    ###    ########       ##     #######     ##    ##     ## ########
//
//                         ######   ########     ###    ########           ##    ##  ######  ########
//                        ##    ##  ##     ##   ## ##   ##     ##          ###   ## ##    ## ##     ##
//                        ##        ##     ##  ##   ##  ##     ##          ####  ## ##       ##     ##
//                        ##   #### ########  ##     ## ##     ##          ## ## ##  ######  ##     ##
//                        ##    ##  ##   ##   ######### ##     ##          ##  ####       ## ##     ##
//                        ##    ##  ##    ##  ##     ## ##     ##          ##   ### ##    ## ##     ##
//                         ######   ##     ## ##     ## ########           ##    ##  ######  ########
void SparseRegressionProblem :: generateCuts_total_gradient_adaptive_nsd( IloEnv &env, IloConstraintArray &cuts2BeAdded, DataModel &data, TuneParameters &tune, IloBoolVarArray &z, IloNumVar &u, IloNumVar &l, IloNumArray valz, IloNum valu, IloNum vall, long nnode, vector<double>  &lastValz)
{
    // backup as copy2
#if test_adap_nsd_cut
    cout << " SparseRegressionProblem :: generateCuts_total_gradient_adaptive_nsd tag (fbw4j7r56hns4e5rtdyc4awe) " << endl;
#endif
    int n = data.getCols();
    
    
    
    double fractional_part = 0;
    for (int i = 0; i < n; i++){
        fractional_part += min( abs(valz[i]), abs(valz[i]-1) );
    }
    if (fractional_part < 1e-4 * n)
        return;
    
    double distFromLastValz = 0;
    for (int i = 0; i < n ; i++){
        ;
        distFromLastValz += pow(lastValz[i] - valz[i], 2);
    }
    cout << " distFromLastValz = " << distFromLastValz << endl;
    if (distFromLastValz < .01*n)
        return;
    
    
    
    double epsilon = .05;
    IloNumArray centralizedValz(env,n);
    for (int i = 0; i < n; i++){
        // find the scale from centralized z.
        // find the gradient from the read z
        if (valz[i] < epsilon )
            centralizedValz[i] = epsilon;
        else if (valz[i] > 1.0-epsilon)
            centralizedValz[i] = 1.0-epsilon;
        else
            centralizedValz[i] = valz[i];
    }
#if test_adap_nsd_cut
    cout << " SparseRegressionProblem :: generateCuts_total_gradient_adaptive_nsd tag (aw4bfyigv4erd) " << endl;
#endif

    vector<double> oneMinusXhatOverXhat(n);
    for (int i = 0; i < n; i++)
        oneMinusXhatOverXhat[i] = (1-centralizedValz[i])/centralizedValz[i];
    vector<double> sqrtOneMinusXhatOverXhat(n);
    for (int i = 0; i < n; i++)
        sqrtOneMinusXhatOverXhat[i] = sqrt(oneMinusXhatOverXhat[i]);

    vector<vector<double>> XinvLX( data.KUK(sqrtOneMinusXhatOverXhat, data.OriginalLinverse) );
    vector<vector<double>> XinvUX( data.KUK(sqrtOneMinusXhatOverXhat, data.OriginalUinverse) );
    
#if test_adap_nsd_cut
    cout << " SparseRegressionProblem :: generateCuts_total_gradient_adaptive_nsd tag (f34nrgwv4b7na4w) " << endl;
#endif

    vector<double> J(n);
    vector<double> K(n);
    for (int i = 0; i < n; i++){
        J[i] = sqrt(1.0/data.OriginalL[i][i]);
        K[i] = sqrt(1.0/data.OriginalU[i][i]);
    }
#if test_adap_nsd_cut || test_adap_nsd_cut_l2
    cout << " tag hrdgerxferbghdrt6nuft ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! " << endl;
    cout << " K = "  << endl;
    for (int i = 0; i < n; i++) cout << K[i] << "\t"; cout << endl;
    cout << " J = "  << endl;
    for (int i = 0; i < n; i++) cout << J[i] << "\t"; cout << endl;
#endif
    

    // fit condition 2
    {
        double maxEigOfKUKplusJLJ = data.eig(data.KUKplusJLJ(K,J,data.OriginalU,data.OriginalL))[n-1];
        double scalerForInitialJK = sqrt(2.0/ maxEigOfKUKplusJLJ);
        for (int i = 0; i < n; i++){
            J[i] *= scalerForInitialJK ;
            K[i] *= scalerForInitialJK ;
        };
    }
#if test_adap_nsd_cut || test_adap_nsd_cut_l2
    cout << " tag xtgnhr6dnse4vges ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! , after fit condition 2" << endl;
    cout << " K = "  << endl;
    for (int i = 0; i < n; i++) cout << K[i] << "\t"; cout << endl;
    cout << " J = "  << endl;
    for (int i = 0; i < n; i++) cout << J[i] << "\t"; cout << endl;
#endif
    
    // fit condition 1
    {
        if (data.eig(data.KUKminusJLJ(K, J, data.OriginalU, data.OriginalL))[n-1] > 1e-15){
            double scalerOfK = data.minEigAsqrtinvBAsqrtinv( data.KUK(K,data.OriginalU),data.KUK(J,data.OriginalL) );
            for (int i = 0; i < n; i++)
                K[i] = K[i]*scalerOfK;
        }
    }
#if test_adap_nsd_cut || test_adap_nsd_cut_l2
    cout << " tag fzsefzresvge54bhr6d ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! , after fit condition 21" << endl;
    cout << " K = "  << endl;
    for (int i = 0; i < n; i++) cout << K[i] << "\t"; cout << endl;
    cout << " J = "  << endl;
    for (int i = 0; i < n; i++) cout << J[i] << "\t"; cout << endl;
#endif
    
    // fit condition 2 again
    {
        double maxEigOfKUKplusJLJ = data.eig(data.KUKplusJLJ(K,J,data.OriginalU,data.OriginalL))[n-1];
        double scalerForInitialJK = sqrt(2.0/ maxEigOfKUKplusJLJ);
        for (int i = 0; i < n; i++){
            J[i] *= scalerForInitialJK ;
            K[i] *= scalerForInitialJK ;
        };
    }
    
    
    
#if test_adap_nsd_cut || test_adap_nsd_cut_l2
    cout << " tag 4f35rvgestf3grthrt ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! , after fit condition 212" << endl;
    cout << " K = "  << endl;
    for (int i = 0; i < n; i++) cout << K[i] << "\t"; cout << endl;
    cout << " J = "  << endl;
    for (int i = 0; i < n; i++) cout << J[i] << "\t"; cout << endl;
#endif

    
    
    
    
#if test_adap_nsd_cut
    cout << " tag w4ehb5dts4e5t ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! " << endl;
    
    cout << " eig = "  << endl;
    for (int i = 0; i < n; i++) cout << data.eig(data.KUKplusJLJ(K,J,data.OriginalU,data.OriginalL))[i] << "\t"; cout << endl;
    cout << " K = "  << endl;
    for (int i = 0; i < n; i++) cout << K[i] << "\t"; cout << endl;
    cout << " J = "  << endl;
    for (int i = 0; i < n; i++) cout << J[i] << "\t"; cout << endl;
#endif

    double deltaTotalObjective = 1e3;
    int maxTotalIterate = 10;
    int iteTotal = 0;
#if test_adap_nsd_cut
    cout << " SparseRegressionProblem :: generateCuts_total_gradient_adaptive_nsd tag (fw3a4v5r65e4awesfx) " << endl;
#endif
    
    
    double totalObjectiveOld;
    {
        vector<double> Ksquare(K.size());
        for (int i = 0; i < n; i++)
            Ksquare[i] = K[i]*K[i];
        vector<double> Jsquare(J.size());
        for (int i = 0; i < n; i++)
            Jsquare[i] = J[i]*J[i];
        totalObjectiveOld = data.logDet_AplusX(Jsquare, XinvLX)-data.logDet_AplusX(Ksquare, XinvUX);
        for(int i =0; i < n; i++){
            totalObjectiveOld += centralizedValz[i]*log(Ksquare[i]) - centralizedValz[i]*log(Jsquare[i]);
        }
    }
#if test_adap_nsd_cut || test_adap_nsd_cut_l2
    cout << " tag f34fawfaw4 ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! " << endl;
    cout << " totalObjectiveOld = " << totalObjectiveOld << endl;
#endif
    
    
    while ( (iteTotal++ < maxTotalIterate) && ( deltaTotalObjective > 1e-1 ) ) {
#if test_adap_nsd_cut
        cout << " tag a3owfj8a4890jfa40 " << endl;
        cout << " iteTotal = " << iteTotal << endl;
        cout << " maxTotalIterate = " << maxTotalIterate << endl;
        cout << " deltaTotalObjective = " << deltaTotalObjective << endl;
#endif
        // move K
        double deltaKObjective = 1e3;
        int iteK = 0;
        int maxKIterate = 10;
        double KstepLength = 5/.618;
#if test_adap_nsd_cut || test_adap_nsd_cut_l2
        cout << " tag aw4v5ge6hrb54dr5  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * " << endl;
        cout << " iteK = " << iteK << endl;
        cout << " maxKIterate = " << maxKIterate << endl;
        cout << " deltaKObjective = " << deltaKObjective << endl;
        cout << " KstepLength = " << KstepLength << endl;
        
#endif

        while ( (iteK++ < maxKIterate) && ( deltaKObjective > 1e-1 ) && ( KstepLength > 1e-2 ) ) {
            vector<double> Ksquare(K.size());
            for (int i = 0; i < n; i++)
                Ksquare[i] = K[i]*K[i];
            
            double KobjectiveOld = - data.logDet_AplusX(Ksquare, XinvUX);
            for(int i =0; i < n; i++){
                KobjectiveOld += centralizedValz[i]*log(Ksquare[i]);
            }

            
            
            vector<double> gradKsquare( data.gradientLDYplusM( Ksquare, XinvUX));
#if test_adap_nsd_cut
            cout << " tag 4nfjtn4zvcrd " << endl;
            cout << " gradKsquare = "  << endl;
            for (int i = 0; i < n; i++) cout << gradKsquare[i] << "\t"; cout << endl;
#endif
            for (int i = 0; i < n; i++)
                gradKsquare[i] = -gradKsquare[i] + centralizedValz[i]/Ksquare[i] ;
#if test_adap_nsd_cut || test_adap_nsd_cut_l2
            cout << " tag aw34fcgse5b43sz4cf4 " << endl;
            cout << " gradKsquare = "  << endl;
            for (int i = 0; i < n; i++) cout << gradKsquare[i] << "\t"; cout << endl;
#endif
            KstepLength = 5.0/.618;
            bool flagK = true;
            while (flagK){
                KstepLength *= .618;
                if (KstepLength < 1e-2)
                    flagK = false;
                else {
                    vector<double> Knew(n);
                    for (int i=0; i<n; i++)
                        Knew[i] = sqrt(Ksquare[i]+gradKsquare[i]*KstepLength);
                    double minKnew = *min_element(Knew.begin(), Knew.end());

#if test_adap_nsd_cut || test_adap_nsd_cut_l2
                    cout << " tag f4e5gdrbhaw3faw4v " << endl;
                    cout << " KstepLength = " << KstepLength << endl;
                    cout << " minKnew = " << minKnew << endl;
                    cout << " K = "  << endl;
                    for (int i = 0; i < n; i++) cout << K[i] << "\t"; cout << endl;
                    cout << " Knew = "  << endl;
                    for (int i = 0; i < n; i++) cout << Knew[i] << "\t"; cout << endl;
                    cout << " Ksquare = "  << endl;
                    for (int i = 0; i < n; i++) cout << Ksquare[i] << "\t"; cout << endl;
                    cout << " gradKsquare = "  << endl;
                    for (int i = 0; i < n; i++) cout << gradKsquare[i] << "\t"; cout << endl;
#endif
                    
                    
                    if (minKnew < 1e-2)
                        ;
                    else {
                        // scale Knew if needed
                        // vector < vector < double > > KUKminusJLJ( const vector < double > &K , const vector < double > &J , const vector < vector < double > > &U, const vector < vector < double > > &L ) ;
#if test_adap_nsd_cut
                        cout << " tag fw4v4esves5 " << endl;
                        cout << " eig = "  << endl;
                        for (int i = 0; i < n; i++) cout << data.eig(data.KUKminusJLJ(Knew, J, data.OriginalU, data.OriginalL))[i] << "\t"; cout << endl;
#endif
                        if (data.eig(data.KUKminusJLJ(Knew, J, data.OriginalU, data.OriginalL))[n-1] > 1e-15){
                            double scalerOfKnew = data.minEigAsqrtinvBAsqrtinv( data.KUK(Knew,data.OriginalU),data.KUK(J,data.OriginalL) );
                            for (int i = 0; i < n; i++)
                                Knew[i] = Knew[i]*scalerOfKnew;
#if test_adap_nsd_cut || test_adap_nsd_cut_l2
                            cout << " tag fevxd5xd5hx " << endl;
                            cout << " scalerOfKnew = " << scalerOfKnew  << endl;
                            cout << " eig = "  << endl;
                            for (int i = 0; i < n; i++) cout << data.eig(data.KUKminusJLJ(Knew, J, data.OriginalU, data.OriginalL))[i] << "\t"; cout << endl;
                            cout << " Knew = "  << endl;
                            for (int i = 0; i < n; i++) cout << Knew[i] << "\t"; cout << endl;
#endif
                        }
//                            // !!!!! to be added: scale Knew
//                            ;
                        vector<double> KnewSquare(Knew.size());
                        for (int i = 0; i < n; i++)
                            KnewSquare[i] = Knew[i]*Knew[i];
                        double KobjectiveNew = -data.logDet_AplusX(KnewSquare, XinvUX);
                        for(int i =0; i < n; i++){
                            KobjectiveNew += centralizedValz[i]*log(KnewSquare[i]);
                        }
#if test_adap_nsd_cut || test_adap_nsd_cut_l2

                        cout << " tag 3e5v4cgaew4vae5 " << endl;
                        cout << " KobjectiveNew = " << KobjectiveNew  << endl;
                        cout << " KobjectiveOld = " << KobjectiveOld  << endl;
#endif
                        
                        
                        if ( KobjectiveNew < KobjectiveOld ) {
                            ;
                        } else {
                            flagK = false;
                            deltaKObjective = KobjectiveNew - KobjectiveOld;
#if test_adap_nsd_cut || test_adap_nsd_cut_l2
                            cout << " tag ferbges5nse5ns5en ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! " << endl;
                            cout << " KstepLength = " << KstepLength << endl;
                            cout << " KobjectiveNew = " << KobjectiveNew << endl;
                            cout << " KobjectiveOld = " << KobjectiveOld << endl;
                            cout << " Knew = "  << endl;
                            for (int i = 0; i < n; i++) cout << Knew[i] << "\t"; cout << endl;
                            cout << " K = "  << endl;
                            for (int i = 0; i < n; i++) cout << K[i] << "\t"; cout << endl;
                            cout << " move K ! watch here: %%%%%% %%%%%% %%%%%% %%%%%% %%%%%% %%%%%% %%%%%% %%%%%% %%%%%% %%%%%% %%%%%% %%%%%% %%%%%% %%%%%% %%%%%%  " << endl;

                            {
                                int temp;
                                cin >> temp;
                            }

                            
#endif
                            for (int i=0; i<n; i++)
                                K[i] = Knew[i];
                            KobjectiveOld = KobjectiveNew;


                        }
                    }
                }
            }
        }
#if test_adap_nsd_cut
        cout << " tag faw4ves4ve4, finished choose K ! ! ! ! ! ! " << endl;
        cout << " K = "  << endl;
        for (int i = 0; i < n; i++) cout << K[i] << "\t"; cout << endl;
        cout << " J = "  << endl;
        for (int i = 0; i < n; i++) cout << J[i] << "\t"; cout << endl;
//        {
//            int temp;
//            cin >> temp;
//        }
#endif
        int maxSIterate = 10;
        double deltaSObjective = 1e3;
        int iteS = 0;
        vector<double> S(n,1.0);
        // test as copy 2
        for (int i = 0; i < n; i++) {
            if (valz[i] < .1)
                S[i]= .03;
            else if (valz[i] >.9)
                S[i] = 30;
            else
                S[i] = 1;
        }
        
        vector<double> Jinverse(n);
        vector<double> Kinverse(n);
        for (int i = 0; i < n; i++){
            Jinverse[i] = 1.0/J[i];
            Kinverse[i] = 1.0/K[i];
        }
        vector<vector<double>> JinvXinvLXJinv(data.KUK( Jinverse, XinvLX ));
        vector<vector<double>> KinvXinvUXKinv(data.KUK( Kinverse, XinvUX ));
        double SstepLength = 1/.618;

        
#if test_adap_nsd_cut
        cout << " tag aw3vaw4vaw4va4w ~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
        cout << " JinvXinvLXJinv = " << endl;
        for (int i = 0; i < n; i++){
            for (int j = 0; j < n; j++)
                cout << JinvXinvLXJinv[i][j] << "\t";
            cout << endl;
        }
        cout << endl;
        
        cout << " KinvXinvUXKinv = " << endl;
        for (int i = 0; i < n; i++){
            for (int j = 0; j < n; j++)
                cout << KinvXinvUXKinv[i][j] << "\t";
            cout << endl;
        }
        cout << endl;
        
        cout << " Jinverse = "  << endl;
        for (int i = 0; i < n; i++) cout << Jinverse[i] << "\t"; cout << endl;
        cout << " Kinverse = "  << endl;
        for (int i = 0; i < n; i++) cout << Kinverse[i] << "\t"; cout << endl;

#endif

        
        
        
        while ( (iteS++ < maxSIterate) && ( deltaSObjective > 1e-2 ) ) { // && ( SstepLength > 1e-2 ) ) {
#if test_adap_nsd_cut || test_adap_nsd_cut_l2
            cout << " tag fa4v4eab4a " << endl;
            cout << " iteS = " << iteS << endl;
            cout << " maxSIterate = " << maxSIterate << endl;
            cout << " deltaSObjective = " << deltaSObjective << endl;
            cout << " SstepLength = " << SstepLength << endl;
            cout << " S = "  << endl;
            for (int i = 0; i < n; i++) cout << S[i] << "\t"; cout << endl;
#endif
  
            double SobjectiveOld = data.logDet_AplusX(S, JinvXinvLXJinv) - data.logDet_AplusX(S, KinvXinvUXKinv);
            double SobjectiveNextS = SobjectiveOld;
            vector<double> gradS1( data.gradientLDYplusM( S, JinvXinvLXJinv));
            vector<double> gradS2( data.gradientLDYplusM( S, KinvXinvUXKinv));
            vector<double> gradS(n);
            for (int i = 0; i < n; i++)
                gradS[i] = gradS1[i] - gradS2[i] ;
            
            double normGradS = 0;
            for (int i = 0; i < n; i++)
                normGradS += gradS[i] * gradS[i] ;
            normGradS = sqrt(normGradS);
            
            
#if test_adap_nsd_cut
            cout << " tag faw4ves4va4wvtes5bvhr" << endl;
            cout << " S = "  << endl;
            for (int i = 0; i < n; i++) cout << S[i] << "\t"; cout << endl;
            cout << " gradS1 = "  << endl;
            for (int i = 0; i < n; i++) cout << gradS1[i] << "\t"; cout << endl;
            cout << " gradS2 = "  << endl;
            for (int i = 0; i < n; i++) cout << gradS2[i] << "\t"; cout << endl;
            cout << " gradS = "  << endl;
            for (int i = 0; i < n; i++) cout << gradS[i] << "\t"; cout << endl;
            
//            {
//                int temp;
//                cin >> temp;
//            }
#endif

            SstepLength = 100/.618;
            bool flagS = true;
            
            while (flagS){
                SstepLength *= .618;
#if test_adap_nsd_cut
                cout << " normGradS * SstepLength = " << normGradS * SstepLength << endl;
#endif
                if ( ( normGradS * SstepLength < 1e-2*n) || (SstepLength < 1e-2) )
                    flagS = false;
                else {
                    vector<double> Snew(n);
                    for (int i=0; i<n; i++)
                        Snew[i] = S[i]+gradS[i]*SstepLength;
                    double minSnew = *min_element(Snew.begin(), Snew.end());
                    
#if test_adap_nsd_cut || test_adap_nsd_cut_l2
                    cout << " tag fwa4vse4bts35b4  ~ ~ ~ ~ ~ " << endl;
                    cout << " SstepLength = " << SstepLength << endl;
                    cout << " minSnew = " << minSnew << endl;
                    cout << " S = "  << endl;
                    for (int i = 0; i < n; i++) cout << S[i] << "\t"; cout << endl;
                    cout << " Snew = "  << endl;
                    for (int i = 0; i < n; i++) cout << Snew[i] << "\t"; cout << endl;
                    cout << " gradS = "  << endl;
                    for (int i = 0; i < n; i++) cout << gradS[i] << "\t"; cout << endl;
                    
                    
                    cout << " \t \t\t\t\t\t\t\t\t\t\t\t check Snew and S and GradS !!! \n watch here: %%%%%% %%%%%% %%%%% %%% % % % %%%%%% %%%% tag favw4g5be5 %% %%%%%% %%%%%% %%%%%% %%%%%% %%%%%% %%%%%% %%%%%% %%%%%%  " << endl;

//                    {
//                        int temp;
//                        cin >> temp;
//                    }
#endif
                    
                    if (minSnew < 1e-2)
                        ;
                    else {
                        vector<double> SnewSqrt(n);
                        for (int i=0; i<n; i++)
                            SnewSqrt[i] = sqrt(Snew[i]);
                        double maxEigSKUKJLJS = data.eig( data.KUK(SnewSqrt, data.KUKplusJLJ(K, J, data.OriginalU, data.OriginalL)) )[n-1]  ;
#if test_adap_nsd_cut
                        cout << " tag fawvf4wbg4eb  ~ ~ ~ ~ ~ " << endl;
                        cout << " SstepLength = " << SstepLength << endl;
                        cout << " maxEigSKUKJLJS = " << maxEigSKUKJLJS << endl;
                        
                        vector<vector<double>> tempM(data.KUKplusJLJ(K, J, data.OriginalU, data.OriginalL));
                        cout << " KUKplusJLJ = " << endl;
                        for (int i = 0; i < n; i++){
                            for (int j = 0; j < n; j++)
                                cout << tempM[i][j] << "\t";
                            cout << endl;
                        }
                        
                        vector<vector<double>> tempM2(data.KUK(SnewSqrt, data.KUKplusJLJ(K, J, data.OriginalU, data.OriginalL)));
                        cout << " S KUKplusJLJ S = " << endl;
                        for (int i = 0; i < n; i++){
                            for (int j = 0; j < n; j++)
                                cout << tempM2[i][j] << "\t";
                            cout << endl;
                        }
                        
                        cout << endl;
                        cout << " eig = "  << endl;
                        for (int i = 0; i < n; i++) cout << data.eig( data.KUK(SnewSqrt, data.KUKplusJLJ(K, J, data.OriginalU, data.OriginalL)) )[i] << "\t"; cout << endl;
                        cout << " S = "  << endl;
                        for (int i = 0; i < n; i++) cout << S[i] << "\t"; cout << endl;
                        cout << " Snew = "  << endl;
                        for (int i = 0; i < n; i++) cout << Snew[i] << "\t"; cout << endl;
                        cout << " SnewSqrt = "  << endl;
                        for (int i = 0; i < n; i++) cout << SnewSqrt[i] << "\t"; cout << endl;
#endif

                        for (int i=0; i<n; i++)
                            Snew[i] = Snew[i]*2/maxEigSKUKJLJS;
                        double SobjectiveNew = data.logDet_AplusX(Snew, JinvXinvLXJinv) - data.logDet_AplusX(Snew, KinvXinvUXKinv);
#if test_adap_nsd_cut|| test_adap_nsd_cut_l2
                        cout << " tag fwefvw4antr6jty   ~ ~ ~ ~ ~ " << endl;
                        
                        cout << " centralizedValz = "  << endl;
                        for (int i = 0; i < n; i++) cout << centralizedValz[i] << "\t"; cout << endl;
                        
                        cout << " SobjectiveNew = " << SobjectiveNew << endl;
                        cout << " SobjectiveOld = " << SobjectiveOld << endl;
                        cout << " Snew = "  << endl;
                        for (int i = 0; i < n; i++) cout << Snew[i] << "\t"; cout << endl;
                        cout << " S = "  << endl;
                        for (int i = 0; i < n; i++) cout << S[i] << "\t"; cout << endl;
//                        {
//                            int temp;
//                            cin >> temp;
//                        }
#endif

                        if ( SobjectiveNew < SobjectiveOld ) {
                            ;
                        } else {
                            flagS = false;
#if test_adap_nsd_cut || test_adap_nsd_cut_l2
                            cout << " tag mf7yiftndr5bs ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! " << endl;
                            cout << " SstepLength = " << SstepLength << endl;
                            cout << " SobjectiveOld = " << SobjectiveOld << endl;
                            cout << " SobjectiveNew = " << SobjectiveNew << endl;
                            cout << " S = "  << endl;
                            for (int i = 0; i < n; i++) cout << S[i] << "\t"; cout << endl;
                            cout << " Snew = "  << endl;
                            for (int i = 0; i < n; i++) cout << Snew[i] << "\t"; cout << endl;
                            cout << " gradS = "  << endl;
                            for (int i = 0; i < n; i++) cout << gradS[i] << "\t"; cout << endl;
                            
//                            cout << " iteS = " << iteS << endl;
//                            cout << " maxSIterate = " << maxSIterate << endl;
//                            cout << " deltaSObjective = " << deltaSObjective << endl;
//                            cout << " SstepLength = " << SstepLength << endl;
                            cout << " move S !!! watch here: %%%%%% %%%%%% %%%%% %%% % % % %%%%%% %%%% tag ihkmy7fmft %% %%%%%% %%%%%% %%%%%% %%%%%% %%%%%% %%%%%% %%%%%% %%%%%%  " << endl;
                            {
                                int temp;
                                cin >> temp;
                            }
                            
#endif
                            for (int i=0; i<n; i++)
                                S[i] = Snew[i];
                            SobjectiveNextS = SobjectiveNew;
                            deltaSObjective = SobjectiveNextS - SobjectiveOld;
                            SobjectiveOld = SobjectiveNextS;
                        }
                    }
                }
            }
            
//            deltaSObjective = SobjectiveNew - SobjectiveOld;
//#if test_adap_nsd_cut || test_adap_nsd_cut_l2
//            cout << " tag f43vf343b5wnt7 ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! " << endl;
//            cout << " SobjectiveNew = " << SobjectiveNew << endl;
//            cout << " SobjectiveOld = " << SobjectiveOld << endl;
//            cout << " S = "  << endl;
//            for (int i = 0; i < n; i++) cout << S[i] << "\t"; cout << endl;
//            cout << " iteS = " << iteS << endl;
//            cout << " maxSIterate = " << maxSIterate << endl;
//            cout << " deltaSObjective = " << deltaSObjective << endl;
//            cout << " SstepLength = " << SstepLength << endl;
//            {
//                int temp;
//                cin >> temp;
//            }
//#endif
//            SobjectiveOld = SobjectiveNew;
        }
        for (int i = 0; i < n; i++){
            double sqrtSi = sqrt(S[i]);
            J[i] = J[i] * sqrtSi;
            K[i] = K[i] * sqrtSi;
        }
        
        {
            vector<double> Ksquare(K.size());
            for (int i = 0; i < n; i++)
                Ksquare[i] = K[i]*K[i];
            vector<double> Jsquare(J.size());
            for (int i = 0; i < n; i++)
                Jsquare[i] = J[i]*J[i];
            double totalObjectiveNew = data.logDet_AplusX(Jsquare, XinvLX)-data.logDet_AplusX(Ksquare, XinvUX);
            for(int i =0; i < n; i++){
                totalObjectiveNew += centralizedValz[i]*log(Ksquare[i]) - centralizedValz[i]*log(Jsquare[i]);
            }
            
            
            deltaTotalObjective = totalObjectiveNew - totalObjectiveOld;
            
#if test_adap_nsd_cut || test_adap_nsd_cut_l2
            cout << " tag fawvfw4cwa, finished choose S ! ! ! ! ! ! " << endl;
            cout << " centralizedValz = "  << endl;
            for (int i = 0; i < n; i++) cout << centralizedValz[i] << "\t"; cout << endl;
            cout << " real z = "  << endl;
            for (int i = 0; i < n; i++) cout << valz[i] << "\t"; cout << endl;
            cout << " S = "  << endl;
            for (int i = 0; i < n; i++) cout << S[i] << "\t"; cout << endl;
            cout << " K = "  << endl;
            for (int i = 0; i < n; i++) cout << K[i] << "\t"; cout << endl;
            cout << " J = "  << endl;
            for (int i = 0; i < n; i++) cout << J[i] << "\t"; cout << endl;

            
            cout << " totalObjectiveNew = " << totalObjectiveNew  << endl;
            cout << " totalObjectiveOld = " << totalObjectiveOld  << endl;
            cout << " deltaTotalObjective = " << deltaTotalObjective  << endl;
            cout << " iteTotal = " << iteTotal  << endl;

            
            {
                int temp;
                cin >> temp;
            }
#endif
            totalObjectiveOld = totalObjectiveNew;
#if test_adap_nsd_cut || test_adap_nsd_cut_l2
            cout << " tag fawefwafw443 ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! " << endl;
            cout << " totalObjectiveOld = " << totalObjectiveOld << ", iteTotal = " << iteTotal << endl;
#endif
            
#if test_adap_nsd_cut || test_adap_nsd_cut_l2
            cout << " tag zsevfaw4bdr6nut7  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * " << endl;
            cout << " iteTotal = " << iteTotal << endl;
            cout << " maxTotalIterate = " << maxTotalIterate << endl;
            cout << " deltaTotalObjective = " << deltaTotalObjective << endl;
            cout << " condition:     while ( (iteTotal++ < maxTotalIterate) && ( deltaTotalObjective > 1e-2 ) )  " << endl;
#endif

            
            
        }

        
        

    }
    
    
    
    
    
    
    
#if test_adap_nsd_cut || test_adap_nsd_cut_l2
    cout << " tag f4wbtes5br65n tujtf ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !  ready to add adaptive nsd cut" << endl;
    cout << " centralizedValz = "  << endl;
    for (int i = 0; i < n; i++) cout << centralizedValz[i] << "\t"; cout << endl;
    cout << " real z = "  << endl;
    for (int i = 0; i < n; i++) cout << valz[i] << "\t"; cout << endl;
    cout << " K = "  << endl;
    for (int i = 0; i < n; i++) cout << K[i] << "\t"; cout << endl;
    cout << " J = "  << endl;
    for (int i = 0; i < n; i++) cout << J[i] << "\t"; cout << endl;
    {
        int temp;
        cin >> temp;
    }
#endif
    
    
    // from here, we can add cut according to J and K
    

    IloRange cut;
    IloExpr rhs(env);
    try {
        vector < vector < double > > scaledL = data.scaleBothSide(data.OriginalL, J);;
        vector < vector < double > > scaledU = data.scaleBothSide(data.OriginalU, K);;
        
#if test_adap_nsd_cut && 0
        cout << " \n\n\n\n tag (7jue5bh63w4v4wsger) " << endl ;
        cout << " data.OriginalL = " << endl ;
        for (int i = 0; i < n; i++){
            for (int j = 0; j < n; j++)
                cout << data.OriginalL[i][j] << "\t" ;
            cout << endl;
        }
        cout << " data.OriginalU = " << endl ;
        for (int i = 0; i < n; i++){
            for (int j = 0; j < n; j++)
                cout << data.OriginalU[i][j] << "\t" ;
            cout << endl;
        }
        cout << " scaledL = " << endl ;
        for (int i = 0; i < n; i++){
            for (int j = 0; j < n; j++)
                cout << scaledL[i][j] << "\t" ;
            cout << endl;
        }
        cout << " scaledU = " << endl ;
        for (int i = 0; i < n; i++){
            for (int j = 0; j < n; j++)
                cout << scaledU[i][j] << "\t" ;
            cout << endl;
        }
        cout << " valz =  " ;
        for (int i = 0; i < n; i++)
            cout << valz[i] << "\t" ;
        cout << endl;
        cout << " J =  " ;
        for (int i = 0; i < n; i++)
            cout << J[i] << "\t" ;
        cout << endl;
        cout << " K =  " ;
        for (int i = 0; i < n; i++)
            cout << K[i] << "\t" ;
        cout << endl;
        
        {
            int temp;
            cin >> temp;
        }
#endif
        
        
        
        double lpart = data.logDet_IplusDiagzW( valz, scaledL )  + log(data.Sigma[n][n]) ;
        double upart = data.logDet_IplusDiagzW( valz, scaledU ) ;
        for (int i = 0; i < n; i++){
            lpart -= 2*log(J[i])*valz[i];
            upart -= 2*log(K[i])*valz[i];
        }
        
        double c = lpart - upart;
        
        /* test if current solution violates the cutting planes. if yes, add the cut */
        double criterior2 = (vall - valu - c);
#if test_adap_nsd_cut || test_adap_nsd_cut_l2
        cout << "   * adaptive nsd cut * * adaptive nsd cut * * adaptive nsd cut * * adaptive nsd cut * * adaptive nsd cut * * adaptive nsd cut * * adaptive nsd cut * * adaptive nsd cut * * adaptive nsd cut *  " << endl;
        cout << "c = " << c << endl;
        cout << "lpart = " << lpart << ", upart = " << upart << ", lpart-upart = " << lpart-upart << endl;
        cout << "vall = " << vall << ", valu = " << valu << ", vall-valu = " << vall-valu << endl;
        cout << " ** criterior2 = (vall - valu - c) = " << criterior2 << endl;
        cout << " valz = "  << valz << endl;
#endif
        
        // rhs should be 0 now;
        if ( ( (nnode < 1000) && ( criterior2 < -1e-2 ) )  || ( (nnode > 1000) && ( criterior2 < -1e-1 ) ) ) {
            rhs += lpart - upart;
            vector<double> grad( data.uMinusLGradient( valz, scaledL, scaledU ) );
            for (int i = 0; i < n; i++){
                grad[i] -= 2*log(J[i]);
                grad[i] += 2*log(K[i]);
            }
            
            for (int i = 0; i < n; i++){
                rhs += grad[i] * ( z[i] - valz[i] );
            }
            
            cut = ( l - u - rhs >= 0 );
#if test_adap_nsd_cut || test_adap_nsd_cut_l2
            cout << " cut = "  << cut << endl;
            {
                int temp;
                cin >> temp;
            }
#endif
            
            cuts2BeAdded.add(cut);
            tune.num_gradient_cut_on_LMinusU_adaptive_nsd++;
        }
        rhs.end();
    } catch (...) {
        rhs.end();
        throw;
    };

    
    
    
    
};




//               ###    ########     ###    ########  ######## #### ##     ## ########    ########  #######  ########    ###    ##
//              ## ##   ##     ##   ## ##   ##     ##    ##     ##  ##     ## ##             ##    ##     ##    ##      ## ##   ##
//             ##   ##  ##     ##  ##   ##  ##     ##    ##     ##  ##     ## ##             ##    ##     ##    ##     ##   ##  ##
//            ##     ## ##     ## ##     ## ########     ##     ##  ##     ## ######         ##    ##     ##    ##    ##     ## ##
//            ######### ##     ## ######### ##           ##     ##   ##   ##  ##             ##    ##     ##    ##    ######### ##
//            ##     ## ##     ## ##     ## ##           ##     ##    ## ##   ##             ##    ##     ##    ##    ##     ## ##
//            ##     ## ########  ##     ## ##           ##    ####    ###    ########       ##     #######     ##    ##     ## ########
//
//             ######   ########     ###    ########     ########   ######  ########
//            ##    ##  ##     ##   ## ##   ##     ##    ##     ## ##    ## ##     ##
//            ##        ##     ##  ##   ##  ##     ##    ##     ## ##       ##     ##
//            ##   #### ########  ##     ## ##     ##    ########   ######  ##     ##
//            ##    ##  ##   ##   ######### ##     ##    ##              ## ##     ##
//            ##    ##  ##    ##  ##     ## ##     ##    ##        ##    ## ##     ##
//             ######   ##     ## ##     ## ########     ##         ######  ########
void SparseRegressionProblem :: generateCuts_total_gradient_adaptive_psd( IloEnv &env, IloConstraintArray &cuts2BeAdded, DataModel &data, TuneParameters &tune, IloBoolVarArray &z, IloNumVar &u, IloNumVar &l, IloNumArray valz, IloNum valu, IloNum vall, long nnode, vector<double>  &lastValz)
{
#if test_adap_psd_cut_level5
    cout << " SparseRegressionProblem :: generateCuts_total_gradient_adaptive_psd tag (fawfawfawfawfthtyj8)  " << endl;
    cout << " tag __FILE__ = " << __FILE__ << ", __LINE__ = " << __LINE__ << endl;
#endif
    int n = data.getCols();
    
    
    
    double fractional_part = 0;
    for (int i = 0; i < n; i++){
        fractional_part += min( abs(valz[i]), abs(valz[i]-1) );
    }
//    if (fractional_part < 1e-4 * n)
//        return;
    
    double distFromLastValz = 0;
    for (int i = 0; i < n ; i++){
        ;
        distFromLastValz += pow(lastValz[i] - valz[i], 2);
        
    }
    cout << " distFromLastValz = " << distFromLastValz << endl;
//    if (distFromLastValz < .01*n)
//        return;

    
    
    
    double epsilon = .05;
    IloNumArray centralizedValz(env,n);
    for (int i = 0; i < n; i++){
        // find the scale from centralized z.
        // find the gradient from the read z
        if (valz[i] < epsilon )
            centralizedValz[i] = epsilon;
        else if (valz[i] > 1.0-epsilon)
            centralizedValz[i] = 1.0-epsilon;
        else
            centralizedValz[i] = valz[i];
    }
    
    vector<double> oneMinusXhatOverXhat(n);
    for (int i = 0; i < n; i++)
        oneMinusXhatOverXhat[i] = (1-centralizedValz[i])/centralizedValz[i];
    vector<double> sqrtOneMinusXhatOverXhat(n);
    for (int i = 0; i < n; i++)
        sqrtOneMinusXhatOverXhat[i] = sqrt(oneMinusXhatOverXhat[i]);
    vector< vector< double > > ConstraintLHS = data.OriginalInvOfInvSum;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            ConstraintLHS[i][j] *= sqrtOneMinusXhatOverXhat[i] * sqrtOneMinusXhatOverXhat[j];

    
    
    
    
    vector<double> y(n);
    for (int i = 0; i < n; i++)
        y[i] = oneMinusXhatOverXhat[i];
    vector<double> sqrty(n);
    for (int i = 0; i < n; i++)
        sqrty[i] = sqrt(y[i]);

    
    
    
    vector< vector< double > > C = data.OriginalInvOfInvSum;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            C[i][j] *= sqrtOneMinusXhatOverXhat[i] * sqrtOneMinusXhatOverXhat[j]; // C = X% * 2 inv[ invU + invL] * X%
    
    vector< vector< double > > LHS = C;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            LHS[i][j] = C[i][j] / sqrty[i] / sqrty[j];
    vector<double> e = data.eig(LHS) ;
    for (int i = 0; i < n; i++)
        y[i] *= e[0]/2.0;
    
    double deltaTotalObjective = 1e3;
    
    vector<double> negativey(n);
    for (int i = 0; i < n; i++)
        negativey[i] = -(y[i]);
    double totalObjectiveOld = data.logDet_AplusX( y , data.OriginalL ) - data.logDet_AplusX( y , data.OriginalU ) + log(data.Sigma[n][n]);
    double mu = 1.0;
    totalObjectiveOld += mu * data.logDet_AplusX( negativey, C );
    for (int i = 0; i < n; i++)
        totalObjectiveOld += mu * log(y[i]);

    while (mu > 1e-6) {
#if test_adap_psd_cut_level5
        cout << " SparseRegressionProblem :: generateCuts_total_gradient_adaptive_psd tag (hrh45h3egers) " << endl;
        cout << "mu = " << mu << endl;
        cout << " tag __FILE__ = " << __FILE__ << ", __LINE__ = " << __LINE__ << endl;
        cout << "\n\n\n\n\n\n\n";
#endif
        bool flagFindOptimalinMu = false;

        
        mu /= 10.0;
    }
    
    
    for (; mu > 1e-6; mu/=10.0) {
#if test_adap_psd_cut_level5
        cout << " SparseRegressionProblem :: generateCuts_total_gradient_adaptive_psd tag (8m8tny4bye8mir67b) " << endl;
        cout << "mu = " << mu << endl;
        cout << " tag __FILE__ = " << __FILE__ << ", __LINE__ = " << __LINE__ << endl;
        cout << "\n\n\n\n\n\n\n";
#endif
        bool flagFindOptimalinMu = false;
        while (!flagFindOptimalinMu) {
            
            vector<double> direction ( data.newtonDirectionPSDSeparationProblem(y, data.OriginalL, data.OriginalU, C, mu) );
            
            double normDirection = 0;
            for (int i = 0; i < n; i++)
                normDirection += direction[i]*direction[i];
            normDirection = sqrt(normDirection);
            if (normDirection < max(1e-1*mu, 1e-3)*n) {
                flagFindOptimalinMu = true;
            } else {
                double YstepLength = 1;
                double YBestStepLength = 1;
                bool flag = true;
//                while ( flag && (YstepLength < 1000) ){
//                    vector <double> newy(n,0.0);
//                    YstepLength *= 3;
//                    for (int i = 0; i < n; i++){
//                        newy[i] = y[i] + YstepLength * direction[i];
//                        if ( newy[i] < 1e-5 ) {
//                            flag = false;
//                            break;
//                        }
//                    }
//                    if (flag) {
//                        vector<vector<double>> matrix(n, vector<double>(n));
//                        for (int i = 0; i < n; i++)
//                            for (int j = 0; j < n; j++)
//                                if (i == j)
//                                    matrix[i][i] = C[i][i]-newy[i];
//                                else
//                                    matrix[i][j] = C[i][j];
//                        if ( data.eig(matrix)[0] < 1e-12 )
//                            flag = false; // otherwise this is a feasible solution
//                        
//                    }
//                    if (flag) {
//                        vector<double> negativey(n);
//                        for (int i = 0; i < n; i++)
//                            negativey[i] = - newy[i] ;
//                        double totalObjectiveNew;
//                        totalObjectiveNew = data.logDet_AplusX( newy , data.OriginalL ) - data.logDet_AplusX( newy , data.OriginalU ) + log(data.Sigma[n][n]);
//                        totalObjectiveNew += mu * data.logDet_AplusX( negativey, C );
//                        for (int i = 0; i < n; i++)
//                            totalObjectiveNew += mu * log(newy[i]);
//                        if ( totalObjectiveNew - totalObjectiveOld < 0.01*YstepLength)
//                            flag = false;
//                        else
//                            YBestStepLength = YstepLength;
//                    }
//                };
//                flag = true;
                YstepLength = YBestStepLength / .5;
                while (flag && (YstepLength > 1e-3 ) ){
#if test_adap_psd_cut_level5
                    cout << " SparseRegressionProblem :: generateCuts_total_gradient_adaptive_psd tag (ktujtders) " << endl;
                    cout << " YstepLength = " << YstepLength << endl;
                    cout << " tag __FILE__ = " << __FILE__ << ", __LINE__ = " << __LINE__ << endl;
#endif
                    bool subFlag = true;
                    YstepLength *= .5;
#if test_adap_psd_cut_level5
                    cout << " SparseRegressionProblem :: generateCuts_total_gradient_adaptive_psd tag (l.ykuyftdhues) " << endl;
                    cout << " direction = " << endl;
                    for (int i = 0; i < n; i++){
                        cout << direction[i] << ", " ;
                    }
                    cout << "\n";

#endif
                    
                    vector < double > newy;
#if test_adap_psd_cut_level5
                    cout << " SparseRegressionProblem :: generateCuts_total_gradient_adaptive_psd tag (ktujtders) " << endl;
                    cout << " tag __FILE__ = " << __FILE__ << ", __LINE__ = " << __LINE__ << endl;
#endif
                    newy.resize(n);
#if test_adap_psd_cut_level5
                    cout << " SparseRegressionProblem :: generateCuts_total_gradient_adaptive_psd tag (34t3ty) " << endl;
                    cout << " tag __FILE__ = " << __FILE__ << ", __LINE__ = " << __LINE__ << endl;
#endif
                    for (int i = 0; i < n; i++){
                        newy[i] = y[i] + YstepLength * direction[i];
                    }
                    for (int i = 0; i < n; i++){
                        if ( newy[i] < 1e-5 ) {
                            subFlag = false;
                            break;
                        }
                    }
#if test_adap_psd_cut_level5
                    cout << " SparseRegressionProblem :: generateCuts_total_gradient_adaptive_psd tag (t43t54y56) " << endl;
#endif
                    if (subFlag) {
                        vector<vector<double>> matrix(n,vector<double>(n));
//                        vector<vector<double>> matrix;
//                        matrix.resize(n);
//                       for (int i = 0; i < n; i++) {
//                           matrix[i].resize(n);
//                        }
                        for (int i = 0; i < n; i++)
                            for (int j = 0; j < n; j++)
                                if (i == j)
                                    matrix[i][i] = C[i][i]-newy[i];
                                else
                                    matrix[i][j] = C[i][j];
                        if ( data.eig(matrix)[0] < 1e-12 )
                            subFlag = false; // otherwise this is a feasible solution
                    }
                    if (subFlag) {
//                        vector<double> negativey(n);
                        for (int i = 0; i < n; i++)
                            negativey[i] = - newy[i] ;
                        double totalObjectiveNew;
                        totalObjectiveNew = data.logDet_AplusX( newy , data.OriginalL ) - data.logDet_AplusX( newy , data.OriginalU ) + log(data.Sigma[n][n]);
                        totalObjectiveNew += mu * data.logDet_AplusX( negativey, C );
                        for (int i = 0; i < n; i++)
                            totalObjectiveNew += mu * log(newy[i]);

                        if (1) { // ( totalObjectiveNew - totalObjectiveOld < 0.01 * YstepLength){
                            ;
                        } else {
                            YBestStepLength = YstepLength;
                            totalObjectiveOld = totalObjectiveNew;
                            flag = false;
                        }
                    }
                };
#if test_adap_psd_cut_level5
                cout << " SparseRegressionProblem :: generateCuts_total_gradient_adaptive_psd tag (fesfbrhyn drb aw) " << endl;
                cout << " flag = " << flag << endl;
                cout << " y = " << endl;
                for (int i = 0; i < n; i++){
                    cout << y[i] << ", " ;
                }
                cout << "\n";

                cout << " tag __FILE__ = " << __FILE__ << ", __LINE__ = " << __LINE__ << endl;
#endif
                if (flag) {
                    // not need to update y
                    flagFindOptimalinMu = true;

                } else {
                    for (int i = 0; i < n; i++)
                        y[i] = y[i] + YstepLength * direction[i];
                    
                }
            }
        }
    }
    
    
    
    
    
    for (int i = 0; i < n; i++)
        sqrty[i] = sqrt(y[i]);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            LHS[i][j] = C[i][j] / sqrty[i] / sqrty[j];
    vector<double> e2 = data.eig(LHS) ;
    for (int i = 0; i < n; i++)
        y[i] *= e2[0];
#if test_adap_psd_cut_level5
    cout << " fjaoseifjaoewjfaewj" << endl;
    cout << " y =  " ;
    for (int i = 0; i < n; i++)
        cout << y[i] << "\t" ;
    cout << endl;
    cout << "oneMinusXhatOverXhat =  " ;
    for (int i = 0; i < n; i++)
        cout << oneMinusXhatOverXhat[i] << "\t" ;
    cout << endl;
#endif


    vector<double> scale(n);
    for (int i = 0; i < n; i++) {
        scale[i] = sqrt( oneMinusXhatOverXhat[i] / y[i] );
    }
    
    
    
//        start to add cut !!!
    /*
    IloRange cut;
    IloExpr rhs(env);
    try {
        vector < vector < double > > scaledL = data.scaleBothSide(data.OriginalL, scale);;
        vector < vector < double > > scaledU = data.scaleBothSide(data.OriginalU, scale);;
        
        
        
#if test_adap_psd_cut_level4
        {
            cout << " \n\n\n\n tag (g6njk8tk5sebe4) " << endl ;
            cout << " data.OriginalL = " << endl ;
            for (int i = 0; i < n; i++){
                for (int j = 0; j < n; j++)
                    cout << data.OriginalL[i][j] << "\t" ;
                cout << endl;
            }
            cout << " data.OriginalU = " << endl ;
            for (int i = 0; i < n; i++){
                for (int j = 0; j < n; j++)
                    cout << data.OriginalU[i][j] << "\t" ;
                cout << endl;
            }
            cout << " scaledL = " << endl ;
            for (int i = 0; i < n; i++){
                for (int j = 0; j < n; j++)
                    cout << scaledL[i][j] << "\t" ;
                cout << endl;
            }
            cout << " scaledU = " << endl ;
            for (int i = 0; i < n; i++){
                for (int j = 0; j < n; j++)
                    cout << scaledU[i][j] << "\t" ;
                cout << endl;
            }
            cout << " valz =  " ;
            for (int i = 0; i < n; i++)
                cout << valz[i] << "\t" ;
            cout << endl;
            cout << " scale =  " ;
            for (int i = 0; i < n; i++)
                cout << scale[i] << "\t" ;
            cout << endl;
            
            cout << " ready to add cut (tag 0af9wj0)  !! !! ! ! ! ! ! !  ! !! ! ! ! ~~ ~~~ ~ ~ ~ ~ ~ ~  ! ! ! !  ! !! ! ! ! ~~ ~~~ ~ ~ ~ ~ ~ ~  ! ! ! !  ! !! ! ! ! ~~ ~~~ ~ ~ ~ ~ ~ ~  ! ! ! !  ! !! ! ! ! ~~ ~~~ ~ ~ ~ ~ ~ ! ! ! ! ~~ ~~~ ~ ~ ~ ~ ~ ~  ! ! ! !  ! !! ! ! ! ~~ ~~~ ~ ~ ~ ~ ~ ~  ! ! ! !  ! !! ! ! ! ~~ ~~~ ~ ~ ~ ~ ~ ~  ! ! ! !  ! !! ! ! ! ~~ ~~~ ~ ~ ~ ~ ~ ~  ! ! ! !  ! !! ! ! ! ~~ ~~~ ~ ~ ~ ~ ~ ~  =  \n" ;

            {
                int temp;
                cin >> temp;
                cout << " ! ! ! !  " ;
                cin >> temp;
            }
        }
#endif

        
        
        double lpart = data.logDet_IplusDiagzW( valz, scaledL )  + log(data.Sigma[n][n]) ;
        double upart = data.logDet_IplusDiagzW( valz, scaledU ) ;
        double c = lpart - upart;
        //  test if current solution violates the cutting planes. if yes, add the cut
        double criterior2 = (vall - valu - c);
#if test_adap_psd_cut_level4
        cout << " \n\n\n\n tag (fawefbs45sb) " << endl ;
        cout << " lpart = " << lpart<< ", upart = " << upart << ", c = " << c << endl;
        cout << " criterior2 = " << criterior2 << endl;
        
#endif
 
        // rhs should be 0 now;
        if ( ( (nnode < 1000) && ( criterior2 < -1e-2 ) )  || ( (nnode > 1000) && ( criterior2 < -1e-1 ) ) ) {
            rhs += lpart - upart;
            vector<double> grad( data.uMinusLGradient( valz, scaledL, scaledU ) );
            for (int i = 0; i < n; i++){
                rhs += grad[i] * ( z[i] - valz[i] );
            }
            
            cut = ( l - u - rhs >= 0 );
#if test_adap_psd_cut
            cout << " cut = "  << cut << endl;
#endif
            
            
            
            // will be back ! ! ! ! !  @ @ @ @ @ 暂时拿掉 @ @ @ @
//            cuts2BeAdded.add(cut);
            tune.num_gradient_cut_on_LMinusU_adaptive_psd++;
        }
        
        
//        rhs.end();

    } catch (...) {
#if test_adap_psd_cut_level4
        cout << " \n\n\n\n tag (fsefwaegaer) " << endl ;
        cout << " ready to add cut !! !! ! ! ! ! ! !  ! !! ! ! ! ~~ ~~~ ~ ~ ~ ~ ~ ~  ! ! ! !  ! !! ! ! ! ~~ ~~~ ~ ~ ~ ~ ~ ~  ! ! ! !  ! !! ! ! ! ~~ ~~~ ~ ~ ~ ~ ~ ~  ! ! ! !  ! !! ! ! ! ~~ ~~~ ~ ~ ~ ~ ~ ~  ! ! ! !  ! !! ! ! ! ~~ ~~~ ~ ~ ~ ~ ~ ~  ! ! ! !  ! !! ! ! ! ~~ ~~~ ~ ~ ~ ~ ~ ~  ! ! ! !  ! !! ! ! ! ~~ ~~~ ~ ~ ~ ~ ~ ~  ! ! ! !  ! !! ! ! ! ~~ ~~~ ~ ~ ~ ~ ~ ~  ! ! ! !  ! !! ! ! ! ~~ ~~~ ~ ~ ~ ~ ~ ~  ! ! ! !  ! !! ! ! ! ~~ ~~~ ~ ~ ~ ~ ~ ~  ! ! ! !  ! !! ! ! ! ~~ ~~~ ~ ~ ~ ~ ~ ~  ! ! ! !  ! !! ! ! ! ~~ ~~~ ~ ~ ~ ~ ! ! ! ! ! ! ! ~~ ~~~ ~ ~ ~ ~ ~ ~  ! ! ! !  ! !! ! ! ! ~~ ~~~ ~ ~ ~ ~ ~ ~  ! ! ! !  ! !! ! ! ! ~~ ~~~ ~ ~ ~ ~ ~ ~  ! ! ! !  ! !! ! ! ! ~~ ~~~ ~ ~ ~ ~ ~ ~  ! ! ! !  ! !! ! ! ! ~~ ~~~ ~ ~ ~ ~ ~ ~  =  " ;
        
        {
            int temp;
            cout << " ! ! ! !  " ;
            cin >> temp;
        }
#endif
//        rhs.end();
        throw;
    };
     */

};



//                         ########  #######  ########    ###    ##
//                            ##    ##     ##    ##      ## ##   ##
//                            ##    ##     ##    ##     ##   ##  ##
//                            ##    ##     ##    ##    ##     ## ##
//                            ##    ##     ##    ##    ######### ##
//                            ##    ##     ##    ##    ##     ## ##
//                            ##     #######     ##    ##     ## ########
//
//                         ######   ########     ###    ########           ##    ##  ######  ########
//                        ##    ##  ##     ##   ## ##   ##     ##          ###   ## ##    ## ##     ##
//                        ##        ##     ##  ##   ##  ##     ##          ####  ## ##       ##     ##
//                        ##   #### ########  ##     ## ##     ##          ## ## ##  ######  ##     ##
//                        ##    ##  ##   ##   ######### ##     ##          ##  ####       ## ##     ##
//                        ##    ##  ##    ##  ##     ## ##     ##          ##   ### ##    ## ##     ##
//                         ######   ##     ## ##     ## ########           ##    ##  ######  ########

void SparseRegressionProblem :: generateCuts_total_gradient_nsd( IloEnv &env, IloConstraintArray &cuts2BeAdded, DataModel &data, TuneParameters &tune, IloBoolVarArray &z, IloNumVar &u, IloNumVar &l, IloNumArray valz, IloNum valu, IloNum vall, long nnode, vector<double>  &lastValz)
{
#if title_everywhere || test_nsd_cut
    cout << " SparseRegressionProblem :: generateCuts_total_gradient_nsd tag (w54bw45b4w) " << endl;
#endif
    int n = data.getCols();
//    if ( (nnode > 1000) && (rand()%10 > 0 ) )
//        return;
    double fractional_part = 0;
    for (int i = 0; i < n; i++){
        fractional_part += min( abs(valz[i]), abs(valz[i]-1) );
    }
//    if (fractional_part < 1e-4 * n)
//        return;
    
    double distFromLastValz = 0;
    for (int i = 0; i < n ; i++){
        ;
        distFromLastValz += pow(lastValz[i] - valz[i], 2);
        
    }
    cout << " distFromLastValz = " << distFromLastValz << endl;
    if (distFromLastValz < .01*n)
        return;
//    double epsilon = .0;
//    IloNumArray realValz(env,n);
//    for (int i = 0; i < n; i++){
//        realValz[i] = valz[i];
//        
//        if (valz[i] < epsilon )
//            valz[i] = epsilon;
//        else if (valz[i] > 1.0-epsilon)
//            valz[i] = 1.0-epsilon;
//    }
    
#if 0 && test_nsd_cut
    {
        cout << "tag foqifj309fj32, input?" << endl;
        int temp;
        cin >> temp;
    }
#endif
    IloRange cut;
    IloExpr rhs(env);
    try {
        double lpart = data.logDet_IplusDiagzW( valz, data.LatNegSemiDef ) + log(data.Sigma[n][n]) ;
        double upart = data.logDet_IplusDiagzW( valz, data.UatNegSemiDef ) ;
        for (int i = 0; i < n; i++)
            lpart -= 2 * log(data.scalerLatNSD[i])*valz[i];
        for (int i = 0; i < n; i++)
            upart -= 2 * log(data.scalerUatNSD[i])*valz[i];
        
        double c = lpart - upart;
        

        
//        for (int i = 0; i < n; i++){
//            c += grad[i] * (realValz[i] - valz[i]);
//        }
        /* test if current solution violates the cutting planes. if yes, add the cut */
        double criterior = ((vall - valu - c)/max(vall, valu));
        double criterior2 = (vall - valu - c);
        
#if test_nsd_cut
        {
            cout << "   * nsd cut * nsd cut * nsd cut * nsd cut * nsd cut * nsd cut * nsd cut * nsd cut * nsd cut * nsd cut * nsd cut * nsd cut * nsd cut * nsd cut * nsd cut * nsd cut * nsd cut * nsd cut * nsd cut * nsd cut * nsd cut * nsd cut *  " << endl;
            cout << "\n\n\n" ;
            cout << " data.logDet_IplusDiagzW( valz, data.LatNegSemiDef ) = " << data.logDet_IplusDiagzW( valz, data.LatNegSemiDef ) ;
            cout << " data.logDet_IplusDiagzW( valz, data.UatNegSemiDef ) = " << data.logDet_IplusDiagzW( valz, data.UatNegSemiDef ) ;
            cout << "\n" ;
            cout << " log(data.Sigma[n][n]) = " << log(data.Sigma[n][n]) << endl;
            cout << " data.scalerLatNSD = \n";
            for (int i = 0; i < n; i++)
                cout << "i = " << i << "," << " data.scalerLatNSD[i] = " << data.scalerLatNSD[i] << "   \n   ";
            cout << " data.scalerUatNSD = \n";
            for (int i = 0; i < n; i++)
                cout << "i = " << i << "," << " data.scalerUatNSD[i] = " << data.scalerUatNSD[i] << "   \n   ";

            for (int i = 0; i < n; i++)
                cout << "i = " << i << "," << " log(data.scalerLatNSD[i]) = " << log(data.scalerLatNSD[i]) << ", valz[i] = " << valz[i] << "   \n   ";
            for (int i = 0; i < n; i++)
                cout << "i = " << i << "," << " log(data.scalerUatNSD[i]) = " << log(data.scalerUatNSD[i]) << ", valz[i] = " << valz[i] << "   \n   ";

            
            cout << "c = " << c << endl;
            cout << "lpart = " << lpart << ", upart = " << upart << ", lpart-upart = " << lpart-upart << endl;
            cout << "vall = " << vall << ", valu = " << valu << ", vall-valu = " << vall-valu << endl;
            cout << "\n\n\n" ;

        }
#endif

#if test_nsd_cut_level2
        cout << "c = " << c << endl;
        cout << "lpart = " << lpart << ", upart = " << upart << ", lpart-upart = " << lpart-upart << endl;
        cout << "vall = " << vall << ", valu = " << valu << ", vall-valu = " << vall-valu << endl;
        cout << "\n\n\n" ;

        cout << " ** criterior = ((vall - valu - c)/max(vall, valu)) = " << criterior << endl;
        cout << " ** criterior2 = (vall - valu - c) = " << criterior2 << endl;
        cout << " valz = "  << valz << endl;
        {
            cout << "tag q2fq34f,  * nsd cut * * nsd cut * * nsd cut * * nsd cut * * nsd cut *  input?" << endl;
            int temp;
            cin >> temp;
        }
        
#endif
        //        if ( criterior2 < - .01 ) {
        if ( ( (nnode < 1000) && ( criterior2 < -1e-2 ) )  || ( (nnode > 1000) && ( criterior2 < -1e-1 ) ) ) {
            vector<double> grad( data.uMinusLGradient( valz, data.LatNegSemiDef, data.UatNegSemiDef ) );
            
//            rhs += lpart - upart;
            rhs += c;
            for (int i = 0; i < n; i++)
                rhs += grad[i] * ( z[i] - valz[i] );
            
            // the linear term
            for (int i = 0; i < n; i++)
                rhs += 2*(log(data.scalerUatNSD[i]) - log(data.scalerLatNSD[i]) ) * ( z[i] - valz[i] );

            cut = ( l - u - rhs >= 0 );
#if test_nsd_cut
            {
                cout << "cut = " << cut << endl;
                cout << "tag fwf34fq, input?" << endl;
                int temp;
                cin >> temp;
            }
#endif

#if pause_total_grad2
            cout << "total cut = " << cut << endl;
            cout << "grad = ["  ;
            for (int i = 0; i < n; i++){
                cout << grad[i] << "\t" ;
            }
            cout << endl;
            
            cout << " valz = "  << valz << endl;
            cout << " vall = " << vall << ", valu = " << valu << endl ;
            
            cout << " lpart = " << lpart;
            cout << " , upart = " << upart;
            cout << " , data.Sigma[n][n] = " << data.Sigma[n][n];
            cout << " , log(data.Sigma[n][n]) = " << log(data.Sigma[n][n]) ;
            cout << " , (vall - valu - c + 0.01 < 0 ? ) = " << vall - valu - c  ;
            cout << " , c = " << c << endl;
            cout << "  ((vall - valu - c)/max(vall, valu))  = " <<  ((vall - valu - c)/max(vall, valu))  << endl;
            
#endif
            
            
            //        if ( criterior < -.0001 ) {
            cuts2BeAdded.add(cut);
            tune.num_gradient_upper_cut_on_LMinusU_nsd++;
#if pause_total_grad2
            cout << " valz = " << valz << endl;
            cout << " vall = " << vall << ", valu = " << valu << ", vall - valu - c = " << vall - valu - c << ", ((vall - valu - c)/max(vall, valu)) = " << ((vall - valu - c)/max(vall, valu)) ;
            cout << ", criterior = " << criterior << ", tf = " << ( criterior < -.0001 ) ;
            cout << endl;
            cout << "current total grad cut = " << cut << endl;
            
            cout << " -- -- each item = " << endl;
            for (int i = 0; i < n; i++)
                cout << "[" << i << "] = " <<  grad[i] * ( z[i] - valz[i] ) << "\t\t";
            cout << endl;
            
#endif
#if pause_total_grad3
            
            cout << " valz = "  << valz << endl;
            cout << " fractional_part of valz = " << fractional_part << endl;
            cout << " vall = " << vall << ", valu = " << valu << ", vall - valu - c = " << vall - valu - c << endl;
            cout << "grad = ["  ;
            for (int i = 0; i < n; i++){
                cout << grad[i] << "\t" ;
            }
            cout << "]" << endl;
            
            
            cout << "current total grad cut = " << cut << endl;
            cout << "\t\t\t\t\t\t\t #total gradient = " << tune.num_gradient_upper_cut_on_LMinusU_nsd;
            cout << ", #total lazy = " << tune.num_lazy_cut;
            cout << ", #total user = " << tune.num_user_cut;
            cout << ", #total sbm_lo_facet = " << tune.num_sbm_lower_cut_on_g;
            cout << ", #total sbm_upp = " << tune.num_sbm_upper_cut_on_f;
            cout  << endl;
            
            if ( 1 ){
                int temp;
                cin >> temp;
            }
            
            
#endif
            
            
#if pause_total_grad4
            
            cout << " valz = "  << valz << endl;
            cout << "current total grad cut = " << cut << endl;
            cout << "\t\t\t\t\t\t\t #total gradient = " << tune.num_gradient_upper_cut_on_LMinusU_nsd;
            cout << ", #total lazy = " << tune.num_lazy_cut;
            cout << ", #total user = " << tune.num_user_cut;
            cout << ", #total sbm_lo_facet = " << tune.num_sbm_lower_cut_on_g;
            cout << ", #total sbm_upp = " << tune.num_sbm_upper_cut_on_f;
            cout  << endl;
            if ( 1 ){
                int temp;
                cin >> temp;
            }
#endif
        }
        
        
#if pauseUpperGrad
        cout << " * const = " << c << endl;
        cout << " * grad = [ ";
        for (int i = 0; i < n; i++)
            cout << grad[i] <<", " ;
        cout << " ]\n";
        
        cout << " * valz = ";
        cout << valz ;
        cout << "\n";
        
        cout << " * cut = { "<< cut << "} "  << endl;
        cout <<" - - (upper gradient) {tag f34fq34f} - - -  "<< endl;
        if ( 1 ){
            int temp;
            cin >> temp;
        }
#endif
        
        
        
#if debugMode_display_cut || debugMode_display_all_cuts
        cout << "\n %% upper gradient cut (fq34fq34fq3) : " << cut << endl;
        //        cout << " @@ currently, valf = " << valf << ", rhs = " << rhs << endl << endl;
#endif
        
#if debugMode_display_callback_sol
        cout << " valz = " << valz;
        cout << ", valf = " << valf << endl;
#endif
        
#if debugMode_display_pause
        {
            cout << " - have a look on the cut above ( ytebe57b )" << endl ;
            int temp;
            cin >> temp;
        }
#endif
        
        rhs.end();
        
    } catch (...) {
        rhs.end();
        
        throw;
    }

};




//                         ########  #######  ########    ###    ##
//                            ##    ##     ##    ##      ## ##   ##
//                            ##    ##     ##    ##     ##   ##  ##
//                            ##    ##     ##    ##    ##     ## ##
//                            ##    ##     ##    ##    ######### ##
//                            ##    ##     ##    ##    ##     ## ##
//                            ##     #######     ##    ##     ## ########
//
//                         ######   ########     ###    ########           ########   ######  ########
//                        ##    ##  ##     ##   ## ##   ##     ##          ##     ## ##    ## ##     ##
//                        ##        ##     ##  ##   ##  ##     ##          ##     ## ##       ##     ##
//                        ##   #### ########  ##     ## ##     ##          ########   ######  ##     ##
//                        ##    ##  ##   ##   ######### ##     ##          ##              ## ##     ##
//                        ##    ##  ##    ##  ##     ## ##     ##          ##        ##    ## ##     ##
//                         ######   ##     ## ##     ## ########           ##         ######  ########

void SparseRegressionProblem :: generateCuts_total_gradient_psd( IloEnv &env, IloConstraintArray &cuts2BeAdded, DataModel &data, TuneParameters &tune, IloBoolVarArray &z, IloNumVar &u, IloNumVar &l, IloNumArray valz, IloNum valu, IloNum vall, long nnode, vector<double>  &lastValz)
{
#if title_everywhere
    cout << " SparseRegressionProblem :: generateCuts_total_gradient_psd tag (g34qg34) " << endl;
#endif
    int n = data.getCols();
    
//    if ( (tune.num_gradient_upper_cut_on_LMinusU > n * 3 ) && (nnode == 0) )
//        return;
//    if ( (nnode > 1000) && (rand()%10 > 0 ) )
//        return;
    
    double fractional_part = 0;
    for (int i = 0; i < n; i++){
        fractional_part += min( abs(valz[i]), abs(valz[i]-1) );
    }
//    if (fractional_part < 1e-4 * n)
//        return;
    
    double distFromLastValz = 0;
    for (int i = 0; i < n ; i++){
        ;
        distFromLastValz += pow(lastValz[i] - valz[i], 2);
    }
    cout << " distFromLastValz = " << distFromLastValz << endl;
    if (distFromLastValz < .01*n)
        return;
    
    double epsilon = .0;
    
    IloNumArray realValz(env,n);
    for (int i = 0; i < n; i++){
        realValz[i] = valz[i];
        
        if (valz[i] < epsilon )
            valz[i] = epsilon;
        else if (valz[i] > 1.0-epsilon)
            valz[i] = 1.0-epsilon;
    }

    
    IloRange cut;
    IloExpr rhs(env);
    // u <= logdet(I + W * diag(z)) , W = Sigma - I
    // u <= logdet(I + W * diag(zhat))  + gradient' * (z - zhat)
    try {
        double lpart = data.logDet_IplusDiagzW( valz, data.L )  + log(data.Sigma[n][n]) ;
        double upart = data.logDet_IplusDiagzW( valz, data.U ) ;
        double c = lpart - upart;
        vector<double> grad( data.uMinusLGradient( valz ) );
        for (int i = 0; i < n; i++){
            c += grad[i] * (realValz[i] - valz[i]);
        }
        


        /* test if current solution violates the cutting planes. if yes, add the cut */
        double criterior = ((vall - valu - c)/max(vall, valu));
        double criterior2 = (vall - valu - c);
#if test_psd_cut
        cout << " SparseRegressionProblem :: generateCuts_total_gradient_psd tag (fawfq34fq34wg34q) " << endl;
        cout << " realValz = " << realValz << endl;
        cout << " valz = " << valz << endl;
        cout << "c = " << c << endl;
        cout << "lpart = " << lpart << ", upart = " << upart << ", lpart-upart = " << lpart-upart << endl;
        cout << " ** criterior2 = (vall - valu - c) = " << criterior2 << endl;
        {
            int temp;
            cin >> temp;
        }
#endif

#if pause_total_grad4 || test_nsd_cut
        cout << "   * psd cut * psd cut * psd cut * psd cut * psd cut * psd cut * psd cut * psd cut * psd cut * psd cut * psd cut * psd cut * psd cut * psd cut * psd cut * psd cut * psd cut * psd cut * psd cut * psd cut *  " << endl;
        cout << "c = " << c << endl;
        cout << "lpart = " << lpart << ", upart = " << upart << ", lpart-upart = " << lpart-upart << endl;
        cout << "vall = " << vall << ", valu = " << valu << ", vall-valu = " << vall-valu << endl;
        
        cout << " ** criterior = ((vall - valu - c)/max(vall, valu)) = " << criterior << endl;
        cout << " ** criterior2 = (vall - valu - c) = " << criterior2 << endl;
        cout << " valz = "  << valz << endl;
#if test_nsd_cut
        {
            cout << "tag fq2fq34gq34, * psd cut * * psd cut * * psd cut * * psd cut * * psd cut * * psd cut * input?" << endl;
            int temp;
            cin >> temp;
        }
#endif
        

#endif
//        if ( criterior2 < - .01 ) {
        if ( ( (nnode < 1000) && ( criterior2 < -1e-2 ) )  || ( (nnode > 1000) && ( criterior2 < -1e-1 ) ) ) {
            rhs += lpart - upart;
//            for (int i = 0; i < n; i++){
//                if ( (grad[i]) > 1e5 )
//                    rhs += 1e5 * ( z[i] - valz[i] );
//                else if ( (grad[i]) < -1e5 )
//                    rhs += -1e5 * ( z[i] - valz[i] );
//                else if ( abs(grad[i]) > 1e-6 )
//                    rhs += grad[i] * ( z[i] - valz[i] );
//            }
//            for (int i = 0; i < n; i++){
//                if ( grad[i] > 1e5 ){
//                    rhs += grad[i] * ( z[i] - valz[i] );
//                } else if ( grad[i] > 1e-5 ){
//                    rhs += grad[i] * ( z[i] - valz[i] ); // check
//                } else if ( grad[i] > 0 ){
//                    rhs += 1e-5 * ( z[i] - valz[i] );
//                } else if ( grad[i] > -1e-5 ){
////                    rhs += grad[i] * ( z[i] - valz[i] );
//                    ;
//                } else if ( grad[i] > -1e5 ){
//                    rhs += grad[i] * ( z[i] - valz[i] ); // check
//                } else {
//                    rhs += 1e5 * ( z[i] - valz[i] );
//                }
//            }
//            for (int i = 0; i < n; i++){
//                if ( grad[i] > 1e1 ){
//                    rhs += 1e1 * ( z[i] - valz[i] );
//                } else if ( grad[i] > 1e-1 ){
//                    rhs += grad[i] * ( z[i] - valz[i] ); // check
//                } else if ( grad[i] > 0 ){
//                    rhs += 1e-1 * ( z[i] - valz[i] );
//                } else if ( grad[i] > -1e-1 ){
//                    //                    rhs += grad[i] * ( z[i] - valz[i] );
//                    ;
//                } else if ( grad[i] > -1e1 ){
//                    rhs += grad[i] * ( z[i] - valz[i] ); // check
//                } else {
//                    rhs += -1e1 * ( z[i] - valz[i] );
//                }
//            }
            for (int i = 0; i < n; i++){
                rhs += grad[i] * ( z[i] - valz[i] );
            }
            
            cut = ( l - u - rhs >= 0 );
#if test_nsd_cut
            {
                cout << " cut = " << cut << endl;
                cout << "tag fq2fq34gq34, * psd cut * * psd cut * * psd cut * * psd cut * * psd cut * * psd cut * input?" << endl;
                int temp;
                cin >> temp;
            }
#endif

#if pause_total_grad2
            cout << "total cut = " << cut << endl;
            cout << "grad = ["  ;
            for (int i = 0; i < n; i++){
                cout << grad[i] << "\t" ;
            }
            cout << endl;
            
            cout << " valz = "  << valz << endl;
            cout << " vall = " << vall << ", valu = " << valu << endl ;
            
            cout << " lpart = " << lpart;
            cout << " , upart = " << upart;
            cout << " , data.Sigma[n][n] = " << data.Sigma[n][n];
            cout << " , log(data.Sigma[n][n]) = " << log(data.Sigma[n][n]) ;
            cout << " , (vall - valu - c + 0.01 < 0 ? ) = " << vall - valu - c  ;
            cout << " , c = " << c << endl;
            cout << "  ((vall - valu - c)/max(vall, valu))  = " <<  ((vall - valu - c)/max(vall, valu))  << endl;
            
#endif

            
//        if ( criterior < -.0001 ) {
            cuts2BeAdded.add(cut);
            tune.num_gradient_upper_cut_on_LMinusU_psd++;
#if pause_total_grad2
            cout << " valz = " << valz << endl;
            cout << " vall = " << vall << ", valu = " << valu << ", vall - valu - c = " << vall - valu - c << ", ((vall - valu - c)/max(vall, valu)) = " << ((vall - valu - c)/max(vall, valu)) ;
            cout << ", criterior = " << criterior << ", tf = " << ( criterior < -.0001 ) ;
            cout << endl;
            cout << "current total grad cut = " << cut << endl;
            
            cout << " -- -- each item = " << endl;
            for (int i = 0; i < n; i++)
                cout << "[" << i << "] = " <<  grad[i] * ( z[i] - valz[i] ) << "\t\t";
            cout << endl;
            
#endif
#if pause_total_grad3
            
            cout << " valz = "  << valz << endl;
            cout << " fractional_part of valz = " << fractional_part << endl;
            cout << " vall = " << vall << ", valu = " << valu << ", vall - valu - c = " << vall - valu - c << endl;
            cout << "grad = ["  ;
            for (int i = 0; i < n; i++){
                cout << grad[i] << "\t" ;
            }
            cout << "]" << endl;
            
            
            cout << "current total grad cut = " << cut << endl;
            cout << "\t\t\t\t\t\t\t #total gradient = " << tune.num_gradient_upper_cut_on_LMinusU_psd;
            cout << ", #total lazy = " << tune.num_lazy_cut;
            cout << ", #total user = " << tune.num_user_cut;
            cout << ", #total sbm_lo_facet = " << tune.num_sbm_lower_cut_on_g;
            cout << ", #total sbm_upp = " << tune.num_sbm_upper_cut_on_f;
            cout  << endl;
            
            if ( 1 ){
                int temp;
                cin >> temp;
            }
            
            
#endif
#if pause_total_grad4
            
//            cout << " valz = "  << valz << endl;
//            cout << " fractional_part of valz = " << fractional_part << endl;
//            cout << " vall = " << vall << ", valu = " << valu << ", vall - valu - c = " << vall - valu - c << endl;
//            cout << "grad = ["  ;
//            for (int i = 0; i < n; i++){
//                cout << grad[i] << "\t" ;
//            }
//            cout << "]" << endl;
//            
            cout << " valz = "  << valz << endl;
            cout << "current total grad cut = " << cut << endl;
            cout << "\t\t\t\t\t\t\t #total gradient = " << tune.num_gradient_upper_cut_on_LMinusU_psd;
            cout << ", #total lazy = " << tune.num_lazy_cut;
            cout << ", #total user = " << tune.num_user_cut;
            cout << ", #total sbm_lo_facet = " << tune.num_sbm_lower_cut_on_g;
            cout << ", #total sbm_upp = " << tune.num_sbm_upper_cut_on_f;
            cout  << endl;
            
            if ( 1 ){
                int temp;
                cin >> temp;
            }
            
            
#endif

            
        }
        
        
#if pauseUpperGrad
        cout << " * const = " << c << endl;
        cout << " * grad = [ ";
        for (int i = 0; i < n; i++)
            cout << grad[i] <<", " ;
        cout << " ]\n";
        
        cout << " * valz = ";
        cout << valz ;
        cout << "\n";
        
        cout << " * cut = { "<< cut << "} "  << endl;
        cout <<" - - (upper gradient) {tag f34fq34f} - - -  "<< endl;
        if ( 1 ){
            int temp;
            cin >> temp;
        }
        
#endif
        
        
        
#if debugMode_display_cut || debugMode_display_all_cuts
        cout << "\n %% upper gradient cut (fq34fq34fq3) : " << cut << endl;
        //        cout << " @@ currently, valf = " << valf << ", rhs = " << rhs << endl << endl;
#endif
        
#if debugMode_display_callback_sol
        cout << " valz = " << valz;
        cout << ", valf = " << valf << endl;
#endif
        
#if debugMode_display_pause
        {
            cout << " - have a look on the cut above ( ytebe57b )" << endl ;
            int temp;
            cin >> temp;
        }
#endif
        
        rhs.end();
        
    } catch (...) {
        rhs.end();
        
        throw;
    }

};


//                    ######## ##     ## ########      ######   ########     ###    ########
//                    ##        ##   ##  ##     ##    ##    ##  ##     ##   ## ##   ##     ##
//                    ##         ## ##   ##     ##    ##        ##     ##  ##   ##  ##     ##
//                    ######      ###    ########     ##   #### ########  ##     ## ##     ##
//                    ##         ## ##   ##           ##    ##  ##   ##   ######### ##     ##
//                    ##        ##   ##  ##           ##    ##  ##    ##  ##     ## ##     ##
//                    ######## ##     ## ##            ######   ##     ## ##     ## ########

void SparseRegressionProblem :: generateCuts_exp_gradient( IloEnv &env, IloConstraintArray &cuts2BeAdded, TuneParameters &tune, IloNumVar &mse, IloNum valmse, IloNumVar &u, IloNum valu, IloNumVar &l, IloNum vall)
{
#if title_everywhere
    cout << " SparseRegressionProblem :: generateCuts_exp_gradient tag (4w5gnw4) " << endl;
#endif

    
    double expX0 = exp(vall - valu);
//    cout << "tag(fawefaw), expX0 = " << expX0 << ", vall = " << vall << ", valu = " << valu << ", valmse = " << valmse << ", (valmse - expX0) = " << (valmse - expX0) << endl;
    if (valmse >= expX0 - 1e-6)
        return;
    IloRange cut;
    IloExpr rhs(env);
    
    rhs = l - u - vall + valu + 1;
    rhs *= expX0;
//    cut = (1993.0 * mse - rhs >= 0);
    cut = (mse - rhs >= 0);
    cuts2BeAdded.add(cut);
    tune.num_gradient_exp++ ;
#if 0 || pause_exp_grad1
    cout << "exp cut = " << cut << endl;
    int temp;
    cin >> temp;
#endif

    rhs.end();
    return;
};



//                ##     ## ########  ########  ######## ########
//                ##     ## ##     ## ##     ## ##       ##     ##
//                ##     ## ##     ## ##     ## ##       ##     ##
//                ##     ## ########  ########  ######   ########
//                ##     ## ##        ##        ##       ##   ##
//                ##     ## ##        ##        ##       ##    ##
//                 #######  ##        ##        ######## ##     ##
//
//                     ######   ########     ###    ########
//                    ##    ##  ##     ##   ## ##   ##     ##
//                    ##        ##     ##  ##   ##  ##     ##
//                    ##   #### ########  ##     ## ##     ##
//                    ##    ##  ##   ##   ######### ##     ##
//                    ##    ##  ##    ##  ##     ## ##     ##
//                     ######   ##     ## ##     ## ########

void SparseRegressionProblem :: generateCuts_uppercuts_gradient( IloEnv &env, IloConstraintArray &cuts2BeAdded, DataModel &data, TuneParameters &tune, IloBoolVarArray &z, IloNumVar &f, IloNumArray &valz, IloNum valf)
{
#if title_everywhere
    cout << " SparseRegressionProblem :: generateCuts_uppercuts_gradient tag (e6m,uwf) " << endl;
#endif

    int n = data.getCols();
    
    IloRange cut;
    IloExpr rhs(env);
//    IloNum valrhs;
#if 0
    IloNum test_rhs = 0;
#endif

    
    
    // g <= logdet(I + W * diag(z)) , W = Sigma - I
    // g <= logdet(I + W * diag(zhat))  + gradient' * (z - zhat)

    try {
        vector<double> grad( data.gradient( valz ) );
        double c = data.logDet_IplusDiagzW( valz ) ;
        rhs += c;
        for (int i = 0; i < n; i++)
            rhs += grad[i] * ( z[i] - valz[i] );
        
        cut = ( f - rhs  <= 0 );
        
#if 0
        {
            cout << " f.getUB() = " << f.getUB() ;
            cout << ",  f.getLB() = " << f.getLB() << endl;
            cout << " c = " << c << ", f = " << valf << endl << endl;
            
            int temp;
            cin >> temp;
        }
#endif

        
        
        if (valf > c + 0.01){
            cuts2BeAdded.add(cut);
            tune.num_gradient_upper_cut_on_f++;
        }
        
        
#if pauseUpperGrad
        cout << " * const = " << c << endl;
        cout << " * grad = [ ";
        for (int i = 0; i < n; i++)
            cout << grad[i] <<", " ;
        cout << " ]\n";
        
        cout << " * valz = ";
        cout << valz ;
        cout << "\n";

        cout << " * cut = { "<< cut << "} "  << endl;
        cout <<" - - (upper gradient) {tag f34fq34f} - - -  "<< endl;
        if ( 1 ){
            int temp;
            cin >> temp;
        }

#endif

        
        
#if debugMode_display_cut || debugMode_display_all_cuts
        cout << "\n %% upper gradient cut (fq34fq34fq3) : " << cut << endl;
//        cout << " @@ currently, valf = " << valf << ", rhs = " << rhs << endl << endl;
#endif
        
#if debugMode_display_callback_sol
        cout << " valz = " << valz;
        cout << ", valf = " << valf << endl;
#endif

#if debugMode_display_pause
        {
            cout << " - have a look on the cut above ( ytebe57b )" << endl ;
            int temp;
            cin >> temp;
        }
#endif
        
        rhs.end();
        
    } catch (...) {
        rhs.end();
        
        throw;
    }
    
};



//                ##     ## ########  ########  ######## ########
//                ##     ## ##     ## ##     ## ##       ##     ##
//                ##     ## ##     ## ##     ## ##       ##     ##
//                ##     ## ########  ########  ######   ########
//                ##     ## ##        ##        ##       ##   ##
//                ##     ## ##        ##        ##       ##    ##
//                 #######  ##        ##        ######## ##     ##
//
//                     ######   ########     ###    ########
//                    ##    ##  ##     ##   ## ##   ##     ##
//                    ##        ##     ##  ##   ##  ##     ##
//                    ##   #### ########  ##     ## ##     ##
//                    ##    ##  ##   ##   ######### ##     ##
//                    ##    ##  ##    ##  ##     ## ##     ##
//                     ######   ##     ## ##     ## ########          max eigen 1

void SparseRegressionProblem :: generateCuts_uppercuts_maxEig1_gradient( IloEnv &env, IloConstraintArray &cuts2BeAdded, DataModel &data, TuneParameters &tune, IloBoolVarArray &z, IloNumVar &f, IloNumArray &valz, IloNum valf)
{
#if title_everywhere || test_upgrad
    cout << " SparseRegressionProblem :: generateCuts_uppercuts_gradient max eigen 1 tag (fq34f3434) " << endl;
    {
        int temp;
        cin >> temp;
    }
#endif
    
    int n = data.getCols();
    
    IloRange cut;
    IloExpr rhs(env);
    //    IloNum valrhs;
    double lx0 = 0;
    
    
    // g <= logdet(I + W * diag(z)) , W = Sigma - I
    // g <= logdet(I + W * diag(zhat))  + gradient' * (z - zhat)
    
    
    
    
    
//    double lpart = data.logDet_IplusDiagzW( valz, data.L )  + log(data.Sigma[n][n]) ;
//    double upart = data.logDet_IplusDiagzW( valz, data.U ) ;
//    double c = lpart - upart;
//    vector<double> grad( data.uMinusLGradient( valz ) );
//    for (int i = 0; i < n; i++){
//        c += grad[i] * (realValz[i] - valz[i]);
//    }
//    
//    
//    
//    /* test if current solution violates the cutting planes. if yes, add the cut */
//    double criterior = ((vall - valu - c)/max(vall, valu));
//    double criterior2 = (vall - valu - c);
//    
//    {
//        cout << "tag fq2fq34gq34, * psd cut * * psd cut * * psd cut * * psd cut * * psd cut * * psd cut * input?" << endl;
//        int temp;
//        cin >> temp;
//    }
//#endif
//    
//    
//#endif
//    //        if ( criterior2 < - .01 ) {
//    if ( ( (nnode < 1000) && ( criterior2 < -1e-2 ) )  || ( (nnode > 1000) && ( criterior2 < -1e-1 ) ) ) {
//        rhs += lpart - upart;
//
//        
    
        
        
        
        
    try {
        double c1 = data.logDet_IplusDiagzW( valz, data.UMaxEigen1 );
        double c2 = 0 ;
        for (int i = 0; i < n; i++)
            c2 += log(data.maxEigenU_ieSigma) * valz[i] ;
        
        double c = c1 + c2;
        
#if test_upgrad
        cout << " in a upper gradient maxEig1 cut . (tag f90wajf02q9jf2340j ) " ;
        cout << " * ME1 cut * ME1 cut * ME1 cut * ME1 cut * ME1 cut * ME1 cut * ME1 cut * ME1 cut * ME1 cut * ME1 cut * ME1 cut *  " << endl ;
        cout << " valz = " << valz << endl;
        cout << " c = " << c << ", valf = " << valf << endl;
//        cout << " data.UMaxEigen1 = " << endl;
//        for (int i = 0; i < n; i++){
//            for (int j = 0; j < n; j++) {
//                cout << data.UMaxEigen1[i][j] << "\t" ;
//            }
//            cout << endl;
//        }
        cout << " data.maxEigenU_ieSigma = " << data.maxEigenU_ieSigma << endl;

        
        if ( c < valf ){ // not cut off
            cout << " going to add cut ! ! ! ! ! " << endl;
        };
        {
            int temp;
            cin >> temp;
        }
#endif
        if ( c >= valf ){ // not cut off
#if test_upgrad
            cout << " no need to add cuts " << endl;
#endif

            rhs.end();
            return;
        };

        rhs += c;
        for (int i = 0; i < n; i++) {
            rhs += ( z[i] - valz[i] ) ;
//            lx0 += valz[i];
        }
        rhs *= log(data.maxEigenU_ieSigma);
//        lx0 *= log(data.maxEigenU_ieSigma);
#if test_upgrad
        cout << " rhs by now = " << rhs << endl;
#endif

        vector < double > grad( data.gradient( valz, data.UMaxEigen1 ) );
        for (int i = 0; i < n; i++)
            rhs += grad[i] * ( z[i] - valz[i] );
        
//        for (int i = 0; i < n; i++)
//            rhs += log(data.maxEigenU_ieSigma) * ( z[i] - valz[i] );
        
        cut = ( f - rhs <= 0 );
        
#if test_upgrad
        cout << " grad = " ;
        for (int i = 0 ; i < n; i++)
            cout << grad[i] << "\t" ;
        cout << "\n";
        cout << " rhs by now = " << rhs << endl;
    
        cout << " cut = " << cut << endl;
        {
            cout << " in me1 cut. tag f90q23j0f9q24j043" << endl;
            int temp;
            cin >> temp;
        }
#endif
        
        if ( valf > c + lx0 + 0.01 ){
            cuts2BeAdded.add(cut);
            tune.num_gradient_upper_cut_on_f++;
        }
        
        
#if pauseUpperGrad
        cout << " * const = " << c << endl;
        cout << " * grad = [ ";
        for (int i = 0; i < n; i++)
            cout << grad[i] <<", " ;
        cout << " ]\n";
        
        cout << " * valz = ";
        cout << valz ;
        cout << "\n";
        
        cout << " * cut = { "<< cut << "} "  << endl;
        cout <<" - - (upper gradient) {tag f34fq34f} - - -  "<< endl;
        if ( 1 ){
            int temp;
            cin >> temp;
        }
        
#endif
        
        
        
#if debugMode_display_cut || debugMode_display_all_cuts
        cout << "\n %% upper gradient cut (fq34fq34fq3) : " << cut << endl;
        //        cout << " @@ currently, valf = " << valf << ", rhs = " << rhs << endl << endl;
#endif
        
#if debugMode_display_callback_sol
        cout << " valz = " << valz;
        cout << ", valf = " << valf << endl;
#endif
        
#if debugMode_display_pause
        {
            cout << " - have a look on the cut above ( ytebe57b )" << endl ;
            int temp;
            cin >> temp;
        }
#endif
        
        rhs.end();
        return;
        
    } catch (...) {
        rhs.end();
        
        throw;
    }
    return;
};




void SparseRegressionProblem :: generateCuts_uppercuts_submodular1( IloEnv &env, IloConstraintArray &cuts2BeAdded, DataModel &data, TuneParameters &tune, IloBoolVarArray &z, IloNumVar &f, IloNumArray &valz, IloNum valf, bitset<_CONST_INT_MAX_FACTOR_SIZE> &currentSelection)
{
#if title_everywhere
    cout << " SparseRegressionProblem :: generateCuts_uppercuts_submodular1 tag (gweg34t3) " << endl;
#endif

    
    int n = data.getCols();
    
    IloRange cut;
    IloExpr rhs(env);
    IloNum valrhs = 0;
#if pause_20140604_0744
    IloNum test_rhs = 0;
#endif
#if pause_20140604_0744
//    cout <<" - - {tag h45wg5wg3w} - - - \n[ IloNum test_rhs; ]\n test_rhs = " << test_rhs << endl;
    
#endif
    /*
     //                             ##     ## ########  ########  ######## ########
     //                             ##     ## ##     ## ##     ## ##       ##     ##
     //                             ##     ## ##     ## ##     ## ##       ##     ##
     //                             ##     ## ########  ########  ######   ########
     //                             ##     ## ##        ##        ##       ##   ##
     //                             ##     ## ##        ##        ##       ##    ##
     //                              #######  ##        ##        ######## ##     ##
     //
     //                              ######  ##     ## ########       ##
     //                             ##    ## ##     ##    ##        ####
     //                             ##       ##     ##    ##          ##
     //                             ##       ##     ##    ##          ##
     //                             ##       ##     ##    ##          ##
     //                             ##    ## ##     ##    ##          ##
     //                              ######   #######     ##        ###### // g
     */
    try {
        bitset<_CONST_INT_MAX_FACTOR_SIZE> fullset;
        for (int j = 0; j < n; j ++)
            fullset.set(j);
        
        rhs += data.logDetSubmatrix( currentSelection );
        valrhs += data.logDetSubmatrix( currentSelection );
#if pause_20140604_0744
        test_rhs += data.logDetSubmatrix( currentSelection );
//        cout <<" - - {tag fq34f3q4t} - - - \n[ test_rhs += data.logDetSubmatrix( currentSelection ); ]\n test_rhs = " << test_rhs << endl;

#endif
        for (int j = 0; j < n; j++){
            //                            cout << " i = " << i << " j = " << j <<  endl;
            bitset<_CONST_INT_MAX_FACTOR_SIZE> jset;
            jset.set(j);
            if (currentSelection[j]){
                rhs += data.logDetDelta( jset, fullset ) * ( z[j] - 1.0 ) ;
                valrhs += data.logDetDelta( jset, fullset ) * ( valz[j] - 1.0 );
#if pause_20140604_0744
                test_rhs += data.logDetDelta( jset, fullset ) * ( test_array[j] - 1.0 );
//                cout << " data.logDetDelta( jset, fullset ) = " << data.logDetDelta( jset, fullset ) << endl;
//                cout <<" - - {tag 3w5g3q4gq34} - - - \n[ test_rhs += data.logDetDelta( jset, fullset ) * ( test_array[j] - 1.0 ); ]\n test_rhs = " << test_rhs << endl;
#endif
            }
            else{
                rhs += data.logDetDelta( jset, currentSelection ) * z[j] ;
                valrhs += data.logDetDelta( jset, currentSelection ) * valz[j] ;
#if pause_20140604_0744
                test_rhs += data.logDetDelta( jset, currentSelection ) * test_array[j] ;
//                cout <<" - - {tag he46h4w5hw453} - - - \n[ test_rhs += data.logDetDelta( jset, currentSelection ) * test_array[j] ; ]\n test_rhs = " << test_rhs << endl;

#endif
            }
        }
        cut = ( f - rhs  <= 0 );
#if pause_20140604_0744
        cout <<" - - (upper cut 1) {tag ofq824n} - - - { > 15.94993} test_rhs = " << test_rhs << endl;
        if ( test_rhs < 15.94993 - 0.1 ){
            cout << " - something wrong here ( tag 2f2323123eq4 )" << endl ;
            int temp;
            cin >> temp;
        }
#endif
        tune.cumulative_violated_amount +=  abs( valf - valrhs ) ;
        
        cuts2BeAdded.add(cut);
        tune.num_sbm_upper_cut_on_f++;

#if debugMode_display_cut || debugMode_display_all_cuts
        cout << " %% upper submodular cut 1 (fq34fq43g5) : " << cut << endl;
//        cout << " @@ currently, valf = " << valf << ", valrhs = " << valrhs << endl << endl;
#endif
#if debugMode_display_callback_sol
        cout << " valz = " << valz;
        cout << ", valf = " << valf << endl;
#endif

        
#if debugMode_display_pause
        {
            cout << " - have a look on the cut above ( 34qpm9cf039 )" << endl ;
            int temp;
            cin >> temp;
        }
#endif

        
        tune.any_cut_added = true;
        tune.num_cut_added_this_round++;
        
        rhs.end();
    } catch (...) {
        rhs.end();
        
        throw;
    }

};

void SparseRegressionProblem :: generateCuts_uppercuts_submodular2( IloEnv &env, IloConstraintArray &cuts2BeAdded, DataModel &data, TuneParameters &tune, IloBoolVarArray &z, IloNumVar &f, IloNumArray &valz, IloNum valf, bitset<_CONST_INT_MAX_FACTOR_SIZE> &currentSelection)
{
#if title_everywhere
    cout << " SparseRegressionProblem :: generateCuts_uppercuts_submodular2 tag (3w5ne564) " << endl;
#endif

    
    int n = data.getCols();

    IloRange cut;
    IloExpr rhs(env);
    IloNum  valrhs = 0;
#if pause_20140604_0744
    IloNum test_rhs = 0;
#endif
     //                             ##     ## ########  ########  ######## ########
     //                             ##     ## ##     ## ##     ## ##       ##     ##
     //                             ##     ## ##     ## ##     ## ##       ##     ##
     //                             ##     ## ########  ########  ######   ########
     //                             ##     ## ##        ##        ##       ##   ##
     //                             ##     ## ##        ##        ##       ##    ##
     //                              #######  ##        ##        ######## ##     ##
     //
     //                              ######  ##     ## ########     #######
     //                             ##    ## ##     ##    ##       ##     ##
     //                             ##       ##     ##    ##              ##
     //                             ##       ##     ##    ##        #######
     //                             ##       ##     ##    ##       ##
     //                             ##    ## ##     ##    ##       ##
     //                              ######   #######     ##       #########
    
    try {
        bitset<_CONST_INT_MAX_FACTOR_SIZE> emptyset;
        rhs += data.logDetSubmatrix( currentSelection );
        valrhs += data.logDetSubmatrix( currentSelection );
#if pause_20140604_0744
        test_rhs += data.logDetSubmatrix( currentSelection );
#endif

        for (int j = 0; j < n; j++){
            bitset<_CONST_INT_MAX_FACTOR_SIZE> jset;
            jset.set(j);
            if (currentSelection[j]){
                rhs += data.logDetDelta( jset, currentSelection )* ( z[j] - 1.0 );
                valrhs += data.logDetDelta( jset, currentSelection )* ( valz[j] - 1.0 );
            }
            else{
                rhs += data.logDetDelta( jset, emptyset ) * z[j];
                valrhs += data.logDetDelta( jset, emptyset ) * valz[j];
            }
        }
        cut = ( f - rhs <= 0 );
        tune.cumulative_violated_amount +=  abs(valf - valrhs) ;

        
        cuts2BeAdded.add(cut);
        tune.num_sbm_upper_cut_on_f++;
        
        
        tune.any_cut_added = true;
        tune.num_cut_added_this_round++;
        rhs.end();
        return;
        
        
    } catch (...) {
        rhs.end();
        
        throw;
    }
};





ILOINCUMBENTCALLBACK6(myIncumbent, IloBoolVarArray,z, IloNumVar,f,  IloNumVar,g, IloNumVar,mse, DataModel &,data, TuneParameters &, tune )
{
#if title_everywhere
    cout << " ILOINCUMBENTCALLBACK6 tag (f34fb34b) " << endl;
#endif
    IloEnv env = getEnv();
    IloNumArray valz(env);
    IloNum valf = 0;
    IloNum valg = 0;
    IloNum valmse = 0;
    try {
        double currentIncumbentValue = getIncumbentObjValue();

        if (currentIncumbentValue < tune.incumbentBestValue){
            tune.incumbentBestValue = currentIncumbentValue;
        };
        
        
        // now we may have a integer feasible solution
        valz  = IloNumArray(env);
        valf  = getValue( f );
        valg  = getValue( g );
        valmse  = getValue( mse );
        getValues( valz, z );
#if pause_at_incumbent
        int n = data.getCols();
        cout << " n = " << n << endl;
        cout << " valz = ["  << endl;
        for (int i = 0; i < n; i++) {
            cout << valz[i] << ",";
        }
        cout << "\n]\n";
        cout << " valf = " << valf << endl;
        cout << " valg = " << valg << endl;
        cout << " mse = " << valmse << endl;
        cout << " \n\n\n here is an incumbent solution!  \n\n\n " << endl;
#endif

        
        valz.end();
        
        int temp;
        cin >> temp;
        
    } catch (...) {
        valz.end();

        throw;
    }

    
    
};


//ILOMIPCallback1(myMIPcallback, TuneParameters &, tune )
//{
//    
//};


ILONODECALLBACK1(myNode, TuneParameters &, tune )
{
//#if test_my_node
//    cout << "getNnodes64() = " << getNnodes64()  << endl;
//    cout << "tune.rootBestValue = " << tune.rootBestValue << endl;
//    cout << "getBestObjValue() = " << getIncumbentObjValue() << endl;
//    {
//        int temp;
//        cin >> temp;
//    }
//#endif
//    if (getNnodes64() == 0){
//        tune.rootBestValue = currentIncumbentValue;
//    };
//

//#if test_my_node
//    cout << "getNnodes64() = " << getNnodes64()  << endl;
//    cout << "tune.rootBestValue = " << tune.rootBestValue << endl;
//    cout << "getBestObjValue() = " << getBestObjValue() << endl;
//    {
//        int temp;
//        cin >> temp;
//    }
//#endif
    if (getNnodes64() == 1){
        tune.rootBestValue = getBestObjValue();
    };
    if (getNnodes64() > tune.limitNumber_Nodes){
        cout << "tune.limitNumber_Nodes = " << tune.limitNumber_Nodes << endl;
        abort();
    };
//        tune.rootBestValue = getBestObjValue();

};





//                    ##          ###    ######## ##    ##        ######  ##     ## ########
//                    ##         ## ##        ##   ##  ##        ##    ## ##     ##    ##
//                    ##        ##   ##      ##     ####         ##       ##     ##    ##
//                    ##       ##     ##    ##       ##          ##       ##     ##    ##
//                    ##       #########   ##        ##          ##       ##     ##    ##
//                    ##       ##     ##  ##         ##          ##    ## ##     ##    ##
//                    ######## ##     ## ########    ##           ######   #######     ##

ILOLAZYCONSTRAINTCALLBACK7(myLazy, IloBoolVarArray,z, IloNumVar,f, IloNumVar,g, IloNumVar,mse, DataModel &,data, TuneParameters &, tune, vector<double> &, lastValz )
{
#if title_everywhere
    cout << " ILOLAZYCONSTRAINTCALLBACK6 tag (gq34b35b) " << endl;
#endif
    if (getNnodes64() == 0){
        tune.rootBestValue = getBestObjValue();
    };

    
    tune.any_cut_added = false;
    tune.num_cut_added_this_round = 0;
    tune.cumulative_violated_amount = 0;
    
#if debugMode_positionTitle
    cout << " \n * * \n * * lazy cut callback: upper cut 12 + lower facet (tag t24vt523v2)  " << endl;
#endif
    
#if debug_20140911_1pm_findSegmentationFault
    cout << "  tag (tq34bw35nw45) " << endl;
#endif
    
    
    //    tune.how_many_usercut_in_one_incumbent = 4;
    
#if debugMode_display_callback_info
    cout << "  getNnodes64 = " << getNnodes64() << ", getNremainingNodes64 = " << getNremainingNodes64() << endl;
    //    cout << "tune.how_many_usercut_in_one_incumbent = " << tune.how_many_usercut_in_one_incumbent << endl;
#endif
    
    IloEnv env = getEnv();
    IloNumArray valz(env);
    IloNum valf = 0;
    IloNum valg = 0;
    IloNum valmse = 0;
    
    IloConstraintArray cuts2BeAddedOnU(env);
    IloConstraintArray cuts2BeAddedOnL(env);
    IloConstraintArray cuts2BeAddedOnLandU(env);

#if debug_20140911_1pm_findSegmentationFault
    cout << "  tag (fw4b w4) " << endl;
#endif
    

    try {
        // now we may have a integer feasible solution
        valz  = IloNumArray(env);
        valf  = getValue( f );
        valg  = getValue( g );
        valmse  = getValue( mse );
        getValues( valz, z );
        
        int n = data.getCols();
#if debugMode_display_callback_sol
        cout << " valz = " << valz;
        cout << ", valf = " << valf;
        cout << ", valg = " << valg << endl;
#endif
#if debug_20140911_1pm_findSegmentationFault
        cout << "  tag (fq34fq34f) " << endl;
#endif
        

        bitset <_CONST_INT_MAX_FACTOR_SIZE> currentSelection;
        for (int j = 0; j < n; j++)
            if ( valz[j] > 0.01 )
                currentSelection.set(j);
#if debugMode_display_sol // debugMode_display_sol
        cout << "    sol information (tag 89uqnoxj034q): " << endl;
        cout << ", currentSelection = " << currentSelection << endl;
        cout << "    sol information (tag 89uqnoxj034q): " << endl;
        

        cout << "    valf = " << valf ;
        double temp = data.logDetSubmatrix(currentSelection) ;
        cout << ", data.logDetSubmatrix(currentSelection) = " << temp << endl;
        cout << "    valg = " << valg ;
        double temp2 = data.logDetSubmatrix_with_y(currentSelection);
        cout << ", data.logDetSubmatrix_with_y(currentSelection) = " << temp2 << endl;
#endif
#if debug_20140911_1pm_findSegmentationFault
        cout << "  tag (f43qfq4f34) " << endl;
        cout << " | valf = " << valf<< endl;
        cout << " | tune.flag_ifuse_LazyConstraints_upper1 = " << tune.flag_ifuse_LazyConstraints_upper1<< endl;
        cout << " | currentSelection = " << currentSelection<< endl;
        cout << " | data.logDetSubmatrix(currentSelection) = ..." << endl;
        cout << " ~ " << data.logDetSubmatrix(currentSelection)<< endl;
        cout << " | tune.DOUBLE_EPSILON_ADD_LAZY_CUT_THRESHOLD = " << tune.DOUBLE_EPSILON_ADD_LAZY_CUT_THRESHOLD<< endl;
        cout << endl;
#endif
        

        // upper cut 1
        if ( ( tune.cutTypeFlag_Lazy_upper1 )
                && ( valf > data.logDetSubmatrix(currentSelection) + tune.DOUBLE_EPSILON_ADD_LAZY_CUT_THRESHOLD ) )
        {
            cout << " apply generateCuts_uppercuts_submodular1 (lazy) \n";
            SparseRegressionProblem :: generateCuts_uppercuts_submodular1(env, cuts2BeAddedOnU, data, tune, z, f, valz, valf, currentSelection);
        }
#if debug_20140911_1pm_findSegmentationFault
        cout << "  tag (gn45s5e) " << endl;
#endif
        

        // upper cut 2
        if ( (tune.cutTypeFlag_Lazy_upper2)
                && ( valf > data.logDetSubmatrix(currentSelection) + tune.DOUBLE_EPSILON_ADD_LAZY_CUT_THRESHOLD ) )
        {
            cout << " apply generateCuts_uppercuts_submodular2 (lazy) \n";
            SparseRegressionProblem :: generateCuts_uppercuts_submodular2(env, cuts2BeAddedOnU, data, tune, z, f, valz, valf, currentSelection);
        }
#if debug_20140911_1pm_findSegmentationFault
        cout << "  tag (w345bw45) " << endl;
#endif
        
        
        
        
        if (tune.cutTypeFlag_Lazy_uppergradient_eig1) {
            cout << " apply generateCuts_uppercuts_maxEig1_gradient (lazy, tag w4fa46v8aw46v) \n";
            SparseRegressionProblem :: generateCuts_uppercuts_maxEig1_gradient(env, cuts2BeAddedOnU, data, tune, z, f, valz, valf);
        }
        

        
        

        if ( ( tune.cutTypeFlag_Lazy_uppergradient )
            && ( valf > data.logDetSubmatrix(currentSelection) + tune.DOUBLE_EPSILON_ADD_LAZY_CUT_THRESHOLD ) )
        {
            cout << " apply generateCuts_uppercuts_gradient (lazy) \n";
            SparseRegressionProblem :: generateCuts_uppercuts_gradient(env, cuts2BeAddedOnU, data, tune, z, f, valz, valf);
        }


        
        
        
        
        
        // lower facet
        if ( ( valg < data.logDetSubmatrix_with_y(currentSelection) - tune.DOUBLE_EPSILON_ADD_LAZY_CUT_THRESHOLD )
                && (tune.cutTypeFlag_Lazy_lowerfacet) ){
            cout << " apply generateCuts_lowercuts (lazy) \n";
            SparseRegressionProblem :: generateCuts_lowercuts(env, cuts2BeAddedOnL, data, tune, z, g, valz, valg);
        }
        
#if debug_20140911_1pm_findSegmentationFault
        cout << "  tag (fq34bw453) " << endl;
#endif

        // exp grad
        if (1){
            cout << " apply generateCuts_exp_gradient (lazy) \n";
            SparseRegressionProblem :: generateCuts_exp_gradient(env, cuts2BeAddedOnL, tune, mse, valmse, f, valf, g, valg);
            cout << " end generateCuts_exp_gradient (lazy) \n";
        }
        //    void static generateCuts_exp_gradient( IloEnv &env, IloConstraintArray &cuts2BeAdded, TuneParameters &tune, IloNumVar &mse, IloNum valmse, IloNumVar &f, IloNum valf, IloNumVar &g, IloNum valg);

        
        
        if  ( tune.cutTypeFlag_Lazy_NSD ) {
            cout << " apply generateCuts_total_gradient_nsd (lazy, tag fawfeawfaewfaw ) \n";
            SparseRegressionProblem :: generateCuts_total_gradient_nsd(env, cuts2BeAddedOnLandU, data, tune, z, f, g, valz, valf, valg, getNnodes(), lastValz);
            cout << " end generateCuts_total_gradient_nsd (lazy, tag fawfeawfaewfaw ) \n";
        }
        
        if  ( tune.cutTypeFlag_Lazy_PSD ) {
            cout << " apply generateCuts_total_gradient_psd (lazy, tag 2r09j2q0jfq) \n";
            SparseRegressionProblem :: generateCuts_total_gradient_psd(env, cuts2BeAddedOnLandU, data, tune, z, f, g, valz, valf, valg, getNnodes(), lastValz);
            cout << " end generateCuts_total_gradient_psd (lazy, tag 2r09j2q0jfq) \n";
        }

        
        
        

#if debugMode_display_cut_info_summary
        cout << " number of cuts added in this round of LAZY on U = " << cuts2BeAddedOnU.getSize() << endl;
        cout << " number of cuts added in this round of LAZY on L = " << cuts2BeAddedOnL.getSize() << endl;
        cout << " number of cuts added in this round of LAZY on L+U = " << cuts2BeAddedOnLandU.getSize() << endl;
//        for (int i = 0; i < cuts2BeAdded.getSize(); i++){
//            cout << "cut[" << i << "] = " << (cuts2BeAdded[i]) << endl;
//        }
#endif
        
        
        for (int i = 0; i < cuts2BeAddedOnU.getSize(); i++){
            if (n < 5)
                cout << cuts2BeAddedOnU[i] << endl;
            
            add( cuts2BeAddedOnU[i] ).end();
            tune.num_lazy_cut++;
            tune.num_lazy_cut_on_f++;
        }
        
        for (int i = 0; i < cuts2BeAddedOnL.getSize(); i++){
            if (n < 5)
                cout << cuts2BeAddedOnL[i] << endl;

            add( cuts2BeAddedOnL[i] ).end();
            tune.num_lazy_cut++;
            tune.num_lazy_cut_on_g++;
        }
        
        for (int i = 0; i < cuts2BeAddedOnLandU.getSize(); i++){
            if (n < 5)
                cout << cuts2BeAddedOnLandU[i] << endl;

            add( cuts2BeAddedOnLandU[i] ).end();
            tune.num_lazy_cut++;
            tune.num_lazy_cut_on_fandg++;
        }
        
        valz.end();
        cuts2BeAddedOnL.end();
        cuts2BeAddedOnU.end();
        cuts2BeAddedOnLandU.end();
#if 1
        cout << " end of current lazy round. tag (faw4ga4gw45) " << endl;
#endif
    } catch (...) {
#if debug_20140911_1pm_findSegmentationFault
        cout << "  tag (f34bq34b) " << endl;
#endif
       valz.end();
        for (int i = 0; i < cuts2BeAddedOnU.getSize(); i++)
            cuts2BeAddedOnU[i].end();
        for (int i = 0; i < cuts2BeAddedOnL.getSize(); i++)
            cuts2BeAddedOnL[i].end();
        for (int i = 0; i < cuts2BeAddedOnLandU.getSize(); i++)
            cuts2BeAddedOnLandU[i].end();
        
        cuts2BeAddedOnLandU.end();
        cuts2BeAddedOnL.end();
        cuts2BeAddedOnU.end();
        throw;
    }
#if debug_20140911_1pm_findSegmentationFault
    cout << "  tag (f34bw345) " << endl;
#endif

#if debugMode_display_pause
    {
        cout << " - at the end of lazy constraint callback" << endl ;
        cout << " - type something here: ( fq34fq4 )" << endl ;
        int temp;
        cin >> temp;
    }
#endif
    
    
};





//
//                    ##     ##  ######  ######## ########         ######  ##     ## ########
//                    ##     ## ##    ## ##       ##     ##       ##    ## ##     ##    ##
//                    ##     ## ##       ##       ##     ##       ##       ##     ##    ##
//                    ##     ##  ######  ######   ########        ##       ##     ##    ##
//                    ##     ##       ## ##       ##   ##         ##       ##     ##    ##
//                    ##     ## ##    ## ##       ##    ##        ##    ## ##     ##    ##
//                     #######   ######  ######## ##     ##        ######   #######     ##

ILOUSERCUTCALLBACK6(myUser, IloBoolVarArray,z, IloNumVar,f, IloNumVar,g, DataModel &,data, TuneParameters &, tune, vector<double> &, lastValz)
{
//    if (getNnodes64() > 10)
//        return;
#if title_everywhere
    cout << " ILOUSERCUTCALLBACK6 tag (6em64en) " << endl;
#endif
    
    if (getNnodes64() == 0){
        tune.rootBestValue = getBestObjValue();
    };

    
    
    
    tune.any_cut_added = false;
    tune.num_cut_added_this_round = 0;
    tune.cumulative_violated_amount = 0;
    
#if debugMode_positionTitle
    cout << " \n * * \n * * user cut callback: upper gradient + lower facet (tag f234v57henu68r)  " << endl;
#endif

#if debugMode_display_callback_info
    cout << "  getNnodes64 = " << getNnodes64() << ", getNremainingNodes64 = " << getNremainingNodes64() << endl;
#endif
    IloEnv env = getEnv();
    IloNumArray valz(env);
    IloNum valf = 0;
    IloNum valg = 0;
    
    IloConstraintArray cuts2BeAddedOnU(env);
    IloConstraintArray cuts2BeAddedOnL(env);
    IloConstraintArray cuts2BeAddedOnLandU(env);

    
    try {
        // now we may have a integer feasible solution
        valz  = IloNumArray(env);
        valf  = getValue( f );
        valg  = getValue( g );
        getValues( valz, z );
        
        int n = data.getCols();
#if debugMode_display_callback_sol
        cout << " valz = " << valz;
        cout << ", valf = " << valf;
        cout << ", valg = " << valg << endl;
#endif
        
        bitset<_CONST_INT_MAX_FACTOR_SIZE> currentSelection;
//        bitset<_CONST_INT_MAX_FACTOR_SIZE> currentNode;
        for (int j = 0; j < n; j++)
            if ( valz[j] > 0.5 )
                currentSelection.set(j);
#if debugMode_display_sol // debugMode_display_sol
        cout << "    sol information (tag fq4293jfq3400): " << endl;
        cout << ",   currentSelection = " << currentSelection << endl;
        cout << "    valf = " << valf << ", data.logDetSubmatrix(currentSelection) = " << data.logDetSubmatrix(currentSelection) << endl;
        cout << "    valg = " << valg << ", data.logDetSubmatrix_with_y(currentSelection) = " << data.logDetSubmatrix_with_y(currentSelection) << endl;
#endif
        
        
        
#if debugMode_usercut && 0
        {
            cout << " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - " << endl ;
            cout << " - - - - - a new label to test when to add user gradient cut" << endl<< endl<< endl ;
            
            cout << " f.getUB() = " << f.getUB() ;
            cout << ",  f.getLB() = " << f.getLB() << endl;
            
            double c = data.logDet_IplusDiagzW( valz ) ;
            cout << " c = " << c << ", f = " << valf << endl << endl;

            int temp;
            cin >> temp;
        }
#endif
        
        // test
        cout << " valz = " <<  valz << endl;
        //
        if  ( tune.cutTypeFlag_User_adaptive_NSD ) {
            cout << " apply generateCuts_total_gradient_adaptive_nsd (user, tag fq4brste4wafers) \n";
            SparseRegressionProblem :: generateCuts_total_gradient_adaptive_nsd(env, cuts2BeAddedOnLandU, data, tune, z, f, g, valz, valf, valg, getNnodes(), lastValz);
            cout << " end generateCuts_total_gradient_adaptive_nsd (user, tag fa3w4egvrv545dt) \n";
        }
        
        if  ( tune.cutTypeFlag_User_NSD ) {
            cout << " apply generateCuts_total_gradient_nsd (user, tag fawfeawfaewfaw ) \n";
            SparseRegressionProblem :: generateCuts_total_gradient_nsd(env, cuts2BeAddedOnLandU, data, tune, z, f, g, valz, valf, valg, getNnodes(), lastValz);
        }
        
        
        
        if  ( tune.cutTypeFlag_User_adaptive_PSD ) {
            cout << " apply generateCuts_total_gradient_adaptive_psd (user, tag fawrbs5enhdr6n) \n";
            SparseRegressionProblem :: generateCuts_total_gradient_adaptive_psd(env, cuts2BeAddedOnLandU, data, tune, z, f, g, valz, valf, valg, getNnodes(), lastValz);
            cout << " end generateCuts_total_gradient_adaptive_psd (user, tag fzsebfdr6m7) \n";
        }
        
        
        if  ( tune.cutTypeFlag_User_PSD ) {
            cout << " apply generateCuts_total_gradient_psd (user, tag fawefawfawgae) \n";
            SparseRegressionProblem :: generateCuts_total_gradient_psd(env, cuts2BeAddedOnLandU, data, tune, z, f, g, valz, valf, valg, getNnodes(), lastValz);
        }
        
        

//        double c = data.logDet_IplusDiagzW( valz ) ;

        
        if (tune.cutTypeFlag_User_uppergradient_eig1) {
            cout << " apply generateCuts_uppercuts_maxEig1_gradient (user, tag 2039f203j) \n";
            SparseRegressionProblem :: generateCuts_uppercuts_maxEig1_gradient(env, cuts2BeAddedOnU, data, tune, z, f, valz, valf);
        }
        
        
        
        // upper gradient
//        if ( ( tune.cutTypeFlag_User_uppergradient )
//                && ( valf > c + tune.DOUBLE_EPSILON_ADD_USER_CUT_THRESHOLD ) ) {
        if (  tune.cutTypeFlag_User_uppergradient  ) {
            cout << " apply generateCuts_uppercuts_gradient (user, tag 2309r23jjf2) \n";
            SparseRegressionProblem :: generateCuts_uppercuts_gradient(env, cuts2BeAddedOnU, data, tune, z, f, valz, valf);
        }
        
        
        
    
        // lower facet
        if ( (tune.cutTypeFlag_User_lowerfacet) ){
            cout << " apply generateCuts_lowercuts (user) \n";
            SparseRegressionProblem :: generateCuts_lowercuts(env, cuts2BeAddedOnL, data, tune, z, g, valz, valg);
        }
        
        
        if  ( 0 ) {
            cout << " apply generateCuts_partial_r_gradient (user) \n";
            SparseRegressionProblem :: generateCuts_partial_r_gradient(env, cuts2BeAddedOnLandU, data, tune, z, f, g, valz, valf, valg, getNnodes(), lastValz);
        }

        
        
        
#if debugMode_display_cut_info_summary
        cout << " number of cuts added in this round of USER on f = " << cuts2BeAddedOnU.getSize() << endl;
        cout << " number of cuts added in this round of USER on g = " << cuts2BeAddedOnL.getSize() << endl;
        cout << " number of cuts added in this round of USER on f-g = " << cuts2BeAddedOnLandU.getSize() << endl;
#endif
        
        for (int i = 0; i < cuts2BeAddedOnU.getSize(); i++){
#if debugMode_display_cut_info_summary
            cout << " in cuts2BeAddedOnU loop, i = " << i << endl;
#endif
            if (n < 5)
                cout << cuts2BeAddedOnU[i] << endl;

            add( cuts2BeAddedOnU[i] ).end();
            tune.num_user_cut++;
            tune.num_user_cut_on_f++;
        }
        
        for (int i = 0; i < cuts2BeAddedOnL.getSize(); i++){
#if debugMode_display_cut_info_summary
            cout << " in cuts2BeAddedOnL loop, i = " << i << endl;
#endif
            if (n < 5)
                cout << cuts2BeAddedOnL[i] << endl;

            add( cuts2BeAddedOnL[i] ).end();
            tune.num_user_cut++;
            tune.num_user_cut_on_g++;
        }
        for (int i = 0; i < cuts2BeAddedOnLandU.getSize(); i++){
#if debugMode_display_cut_info_summary
            cout << " in cuts2BeAddedOnLandU loop, i = " << i << endl;
#endif
            if (n < 5)
                cout << cuts2BeAddedOnLandU[i] << endl;
            
            add( cuts2BeAddedOnLandU[i] ).end();
            tune.num_user_cut++;
            tune.num_user_cut_on_fandg++;
        }
        
#if debugMode_display_pause || debugMode_usercut
        cout << " valz = " << valz;
        cout << ", valf = " << valf;
        cout << ", valg = " << valg << endl;
#endif
        
        
        
        for (int i = 0; i < n ; i++)
            lastValz[i] = valz[i];

        valz.end();
        cuts2BeAddedOnL.end();
        cuts2BeAddedOnLandU.end();
        cuts2BeAddedOnU.end();
        
    } catch (...) {
        valz.end();
        for (int i = 0; i < cuts2BeAddedOnU.getSize(); i++)
            cuts2BeAddedOnU[i].end();
        for (int i = 0; i < cuts2BeAddedOnL.getSize(); i++)
            cuts2BeAddedOnL[i].end();
        for (int i = 0; i < cuts2BeAddedOnLandU.getSize(); i++)
            cuts2BeAddedOnLandU[i].end();
        
        cuts2BeAddedOnL.end();
        cuts2BeAddedOnLandU.end();
        cuts2BeAddedOnU.end();
        throw;

    }
#if debugMode_display_cut_info_summary
    {
        cout << "  end of user cut round jwaoejfaw9jfaw4" << endl ;
    }
#endif

#if debugMode_display_pause || debugMode_usercut
    {
        cout << " - at the end of user cut callback" << endl ;
        cout << " - type something here: ( 45hn667rhe )" << endl;
        int temp;
        cin >> temp;
    }
#endif

    
};








void SparseRegressionProblem :: readData (const string & filename)
{
#if title_everywhere
    cout << " SparseRegressionProblem :: readData tag (g4be54b) " << endl;
#endif

#if debugMode_positionTitle
    cout << "  \n * * \n * * readData. tag jn7n67jmue6be4 " << endl;
#endif
    
#if debugMode_read
    cout << " start to read " << filename << endl;
#endif
    Qi_IO_input :: readData(filename, data);
    
#if debugMode_read
    cout << " finish reading " << filename << endl;
    cout << " data.m = " << data.getRows() << ", data.n = " << data.getCols() << endl;
    cout << " first 20 rows: " << endl;
    for (int i = 0; i < min(data.getRows(), 20); ++i)
        cout << data.bitData.at(i) << endl;
    cout << " last  20 rows: " << endl;
    for (int i = max(0, data.getRows() - 20)  ; i < data.getRows(); ++i)
        cout << data.bitData.at(i) << endl;
#endif
    
    
    data.flag_if_sigma_initialized = false;
    
#if 1 // debugMode_pause
    cout << " \n\n\n here is the end of readData! input something:  \n\n\n " << endl;
    int tmp;
    cin >> tmp;
#endif
};





/*
//                  ######   #######  ########  ########
//                 ##    ## ##     ## ##     ## ##
//                 ##       ##     ## ##     ## ##
//                 ##       ##     ## ########  ######
//                 ##       ##     ## ##   ##   ##
//                 ##    ## ##     ## ##    ##  ##
//                  ######   #######  ##     ## ########
 */

void SparseRegressionProblem :: core (bool flagCPX)
{



#if debugMode_display_stage || title_everywhere
    cout << " SparseRegressionProblem :: core (f3vg6n78)" << endl;
#endif
    //  test();
    
    if ( !data.flag_if_sigma_initialized ){
#if debugMode_display_stage
        cout << " data.initializeSigma(); (vq35ge46n)" << endl;
#endif

        data.initializeSigma();
#if debugMode_display_stage
        cout << " data.scalizeSigma(); (h6ue6hbw354b)" << endl;
#endif
        data.scalizeSigma();

        data.initializeUandLatNSD(tune);
        data.initializeOriginalInvOfInvSum(tune);
        
#if debugMode_display_stage || 1
        cout << " end of data.initializeUandLatNSD(); (faw4f4fa)" << endl;
        {
            int temp;
            cin >> temp;
        }
#endif

        data.initializerStar();
    }
    
    n = data.getCols();
    m = data.getRows();
#if 1
    cout << " problem info: (tag 23cc45ybr67mty8) " << endl;
    cout << " m = " << m << ", n = " << n  << endl;
#endif

#if debugMode_display_round1
    cout << " here we are at (tag fqj4239fjq34098j) " << endl;
#endif
//    setTuneParameters();
    lastValz.resize(n, -1);

    try {
        
        IloModel        model(env);
        IloCplex        cplex(env);
        IloBoolVarArray z(env, n);                    // x[i][j] = x[j*n+i]
        IloNumVar  f(env, -5000, 5000);    // H( Pa(Xi) )
        IloNumVar  g(env, -5000, 5000);    // H( Pa(Xi) | Xi)
        
        IloNumVar  mse(env, 1e-6, 1e6);
        
//        m = 1993;
        
        
        build_SR_VarNames ( model, z, f, g , mse ) ;
        build_SR_Objective ( model, z, f, g , mse ) ;
        build_SR_UpperConstraints ( model, z, f ) ;
        build_SR_LowerConstraints ( model, z, g ) ;
        build_SR_ULboth_Constraints( model, z, f, g ) ;
        build_SR_Exp_Constraints( model, mse, f, g ) ;
        build_SR_Cardinality_Constraints ( model, z ) ;
        
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        //
        
        //
        //     ######     ###    ##       ##       ########     ###     ######  ##    ##
        //    ##    ##   ## ##   ##       ##       ##     ##   ## ##   ##    ## ##   ##
        //    ##        ##   ##  ##       ##       ##     ##  ##   ##  ##       ##  ##
        //    ##       ##     ## ##       ##       ########  ##     ## ##       #####
        //    ##       ######### ##       ##       ##     ## ######### ##       ##  ##
        //    ##    ## ##     ## ##       ##       ##     ## ##     ## ##    ## ##   ##
        //     ######  ##     ## ######## ######## ########  ##     ##  ######  ##    ##
        
        //        cplex.use(MyBranch(env, var));
        //        cplex.use(MySelect(env));
        //        cplex.use(myIncumbent( env, x, f, g, p, data, tune ));
        
        
        
        //        cplex.use( myUser_upperEnvelope1( env, x, f, data ) );
        //        cplex.use( myLazy_upperEnvelope1( env, x, f, g, data ) );
        //        cplex.use( myUser_lowerEnvelope1( env, x, g, data ) );
        //        cplex.use( myLazy_lowerEnvelope1( env, x, g, data ) );
        
        cplex.use( myLazy( env, z, f, g, mse, data, tune, lastValz ) );
        cplex.use( myIncumbent( env, z, f, g, mse, data, tune ) );
        cplex.use( myUser( env, z, f, g, data, tune, lastValz) );
        cplex.use( myNode( env, tune ));
        //        cplex.use( myUser_lowerEnvelope1( env, x, g, data, tune ) );
        
        //        cplex.use( myHeuristic_Feed( env, x, f, g, p, data, tune ) );
        
        
        
        // ILOUSERCUTCALLBACK3(myUser_lowerEnvelope1, IloBoolVarArray &,x, IloNumVarArray &,g, DataModel &, data )
        
        
        
        
        
        
        //        cplex.use( myIncumbent_all1( env, x, f, g, data ) );
        
        
        
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        
        /*
                 ########     ###    ########     ###    ##     ##
                 ##     ##   ## ##   ##     ##   ## ##   ###   ###
                 ##     ##  ##   ##  ##     ##  ##   ##  #### ####
                 ########  ##     ## ########  ##     ## ## ### ##
                 ##        ######### ##   ##   ######### ##     ##
                 ##        ##     ## ##    ##  ##     ## ##     ##
                 ##        ##     ## ##     ## ##     ## ##     ##
         */
        //        cplex.setParam(IloCplex::MIPSearch, IloCplex::Traditional);
        cplex.setParam(IloCplex::MIPDisplay, 3);
        cplex.setParam(IloCplex::MIPInterval, 1);
        cplex.setParam(IloCplex::NumericalEmphasis, 1);
        cplex.setParam(IloCplex::ParallelMode, 1);
        cplex.setParam(IloCplex::Threads, 4);
        cplex.setParam(IloCplex::MIPEmphasis, 4);
        //        cplex.setParam(IloCplex::EpGap, .1);   // Any number from 0.0 to 1.0; default: 1e-04.

//        cplex.setParam(IloCplex::EpAGap, 1e-12);  // Any nonnegative number; default: 1e-06.
//        cplex.setParam(IloCplex::EpGap, 1e-12);   // Any number from 0.0 to 1.0; default: 1e-04.
//        cplex.setParam(IloCplex::EpLin, 1e-12);   // Any positive value greater than zero; default: 1e-3.
//        cplex.setParam(IloCplex::EpMrk, 0.99999); // Any number from 0.0001 to 0.99999; default: 0.01.
//        cplex.setParam(IloCplex::EpInt, 0 ); // Any number from 0.0 to 0.5; default: 1e-05.
//        cplex.setParam(IloCplex::EpPer, 1e-8); // Any positive number greater than or equal to 1e-8; default: 1e-6.
//        cplex.setParam(IloCplex::EpRHS, 1e-9); // Any number from 1e-9 to 1e-1; default: 1e-06.
//        cplex.setParam(IloCplex::EpRelax, 1e-12); // Any nonnegative value; default: 1e-6.
//        cplex.setParam(IloCplex::EpOpt, 1e-09); // Any number from 1e-9 to 1e-1; default: 1e-06.
//        
//        
//        cplex.setParam(IloCplex::BtTol, 0); // Any number from 1e-9 to 1e-1; default: 1e-06.
        
        
        //        cplex.setParam(IloCplex::SIM, 1e-09); // Any number from 1e-9 to 1e-1; default: 1e-06.
        
        
        //        cplex.setParam(IloCplex::EpMrk, 1e-6);
        //        cplex.setParam(IloCplex::Ep, 1e-6);
        
        //        cplex.setParam(IloCplex::PreInd, 0);
        // cplex.setParam(IloCplex::MIPEmphasis, 1);
        
        


        

        // this is for cplex 12.6
//        cplex.setParam(IloCplex::Param::MIP::Cuts::Gomory, -1);
//        cplex.setParam(IloCplex::Param::MIP::Cuts::FlowCovers, -1);
//        cplex.setParam(IloCplex::Param::MIP::Cuts::MIRCut, -1);
        
        cplex.extract(model);
        cplex.exportModel ("model1.lp");
        
        cplex.solve();
        cplex.writeSolution("./testWriteSol.sol");
        
        output_SR_Solutions ( cplex, z, f, g, mse, tune ) ;
        
    }
    catch (IloException& e) {
        cerr << "Concert exception caught: " << e << endl;
    }
    catch (...) {
        cerr << "Unknown exception caught" << endl;
    }
    env.end();
    
    
    
    //    test();
    
    
    
};





/*
//              #######  ########        ## ########  ######  ######## #### ##     ##
//             ##     ## ##     ##       ## ##       ##    ##    ##     ##  ##     ##
//             ##     ## ##     ##       ## ##       ##          ##     ##  ##     ##
//             ##     ## ########        ## ######   ##          ##     ##  ##     ##
//             ##     ## ##     ## ##    ## ##       ##          ##     ##   ##   ##
//             ##     ## ##     ## ##    ## ##       ##    ##    ##     ##    ## ##
//              #######  ########   ######  ########  ######     ##    ####    ###
 */

void SparseRegressionProblem :: build_SR_Objective (IloModel &model, IloBoolVarArray &z, IloNumVar &f, IloNumVar &g , IloNumVar &mse )
{
#if title_everywhere
    cout << " SparseRegressionProblem :: build_SR_Objective tag (e465me57) " << endl;
#endif


#if debugMode_positionTitle
    cout << "  \n * * \n * * build_SR_Objective. tag f34vh57hn " << endl;
#endif
    IloNumExpr exprObj(env);
    IloNumExpr exprSumz(env);
    for (int i = 0; i < n; i++) {
        exprSumz += z[i];
    }

//    exprObj = g - f;
    
//    exprObj = m * mse + log(m)*10 * exprSumz;
    exprObj = mse;
    
    cout << "the obj is " << endl;
    cout << exprObj << endl;
    {
        cout << " \n\n\n here is the obj!  input something here: \n\n\n " << endl;

        int temppp;
        cin >> temppp;
    }
//
//    for (int i = 0; i < n; i++) {
//        exprObj += z[i] * 1.64975364935;
//    }
    model.add(IloMinimize(env, exprObj));
    env.out() << exprObj << "\n";
    
    
    
    // cardinality constraints
    if (tune.kSparse == 0) {
        int temppp;
        cout << "illegal ! ! ! !  kSparse = " << tune.kSparse << endl;
        cin >> temppp;
        return;
    };
        int k = tune.kSparse;
        IloExpr  rhs(env);
        for ( int j = 0; j < n; j++) {
            rhs += z[j];
        };
        model.add( rhs <= k);
    

    
    
    
};




void SparseRegressionProblem :: build_SR_Exp_Constraints (IloModel &model, IloNumVar &mse, IloNumVar &u, IloNumVar &l)
{
#if title_everywhere
    cout << " SparseRegressionProblem :: build_SR_Exp_Constraints tag (4meh4nm) " << endl;
#endif

    model.add(mse >= l - u + 1);
    
    cout << "Exp_Constraints = " << ( mse >= l - u + 1) << endl;
    cout << " \n\n\n here is the exp constraint!  input something here: \n\n\n " << endl;
    int temp;
    cin >> temp;
};


void SparseRegressionProblem :: build_SR_VarNames (IloModel &model, IloBoolVarArray &z, IloNumVar &f, IloNumVar &g, IloNumVar &mse  )
{
#if title_everywhere
    cout << " SparseRegressionProblem :: build_SR_VarNames tag (fwbw43b) " << endl;
#endif


#if debugMode_positionTitle
    cout << "  \n * * \n * * build_SR_VarNames. tag f3rbm78 " << endl;
#endif
    int n = data.getCols();
    
    for (int i = 0; i < n; i++) {
        string name = "z";
        name += to_string( i+1 );
        z[i].setName(name.c_str());
#if 0
        cout << z[i] << endl;
#endif
    };
    
    f.setName("u");
    g.setName("l");
    mse.setName("mse");
#if 0
    cout << f << endl;
    cout << g << endl;
#endif
};



void SparseRegressionProblem :: build_SR_ULboth_Constraints (IloModel &model, IloBoolVarArray &z, IloNumVar &u, IloNumVar &l )
{
#if title_everywhere
    cout << " SparseRegressionProblem :: build_SR_ULboth_Constraints tag (34b34n ) " << endl;
#endif
    if (!tune.initialConstraintFlag_totalGrad_PSD)
        return;
    
    
    // 0000
    IloExpr rhs0(env);
    {
        IloNumArray zeros(env, n);
        for (int i = 0; i < n; i++){
            zeros[i] = 0;
        }
        vector<double> grad0( data.uMinusLGradient( zeros ) );
        double lpart0 = data.logDet_IplusDiagzW( zeros, data.L )  + log(data.Sigma[n][n]) ;
        double upart0 = data.logDet_IplusDiagzW( zeros, data.U ) ;
        double c0 = lpart0 - upart0;
        rhs0 += c0;
        for (int i = 0; i < n; i++)
            rhs0 += grad0[i] * ( z[i]);
    }
    model.add( l - u - rhs0 >= 0 );

    
    // 1111
    IloExpr rhs1(env);
    {
        IloNumArray ones(env, n);
        for (int i = 0; i < n; i++){
            ones[i] = 1;
        }
        vector<double> grad1( data.uMinusLGradient( ones ) );
        double lpart1 = data.logDet_IplusDiagzW( ones, data.L )  + log(data.Sigma[n][n]) ;
        double upart1 = data.logDet_IplusDiagzW( ones, data.U ) ;
        double c1 = lpart1 - upart1;
        rhs1 += c1;
        
        for (int i = 0; i < n; i++)
            rhs1 += grad1[i] * ( z[i] - 1.0 );
    }
    model.add( l - u - rhs1 >= 0 );
    
    
    IloExpr rhs2(env);
    {
        IloNumArray kOverNs(env, n);
        double temp = ((double)tune.kSparse) / n;
        for (int i = 0; i < n; i++){
            kOverNs[i] = temp;
        }
        vector<double> grad2( data.uMinusLGradient( kOverNs ) );
        double lpart2 = data.logDet_IplusDiagzW( kOverNs, data.L )  + log(data.Sigma[n][n]) ;
        double upart2 = data.logDet_IplusDiagzW( kOverNs, data.U ) ;
        double c2 = lpart2 - upart2;
        rhs2 += c2;
        
        for (int i = 0; i < n; i++)
            rhs2 += grad2[i] * ( z[i] - kOverNs[i] );
    }
    model.add( l - u - rhs2 >= 0 );


    
#if 1 || pause_initial_cons
    cout << " initial total grad =  " << endl;
    cout << l - u - rhs0 << " >= 0 " << endl;
    cout << l - u - rhs1 << " >= 0 " << endl;
    cout << l - u - rhs2 << " >= 0 " << endl;
    {
        cout << " \n\n\n input something here: \n\n\n " << endl;
        int temp;
        cin >> temp;
    }
#endif

    rhs0.end();
    rhs1.end();
    rhs2.end();
};


void SparseRegressionProblem :: build_SR_UpperConstraints (IloModel &model, IloBoolVarArray &z, IloNumVar &f)
{
    
    if (!tune.initialConstraintFlag_upper)
        return;

#if title_everywhere
    cout << " SparseRegressionProblem :: build_SR_UpperConstraints tag (wa4 na4m ) " << endl;
#endif


#if debugMode_positionTitle
    cout << "  \n * * \n * * build_SR_UpperConstraints. (tag fq34f34) " << endl;
#endif

    int n = data.getCols();
    IloNumExpr rhs(env);
    
    // this is the 20140604 version
//    rhs += 1;
//    for (int j = 0; j < n; j++){
//        bitset<_CONST_INT_MAX_FACTOR_SIZE> emptyset;
//        emptyset.set(n);
//        rhs += data.logDetDelta(j, emptyset) * z[j]; //x[j][i]
//    };

    for (int j = 0; j < n; j++){
//        rhs += ( data.Sigma[j][j] - 1 ) * z[j]; //x[j][i]
        rhs += log( data.Sigma[j][j] ) * z[j]; //x[j][i]
    };

   
    
#if debug_20140911_1pm_findSegmentationFault || debugMode
    cout << "  \n * * \n * * build_SR_UpperConstraints. tag fw34v4eg5b6 " << endl;
#endif
    model.add( f - rhs <= 0 );
    
    
    IloNumExpr rhs2(env);
    bitset<_CONST_INT_MAX_FACTOR_SIZE> fullset;
    for (int j = 0; j < n; j++)
        fullset.set(j);
#if debug_20140911_1pm_findSegmentationFault
    cout << " tag f23f23f234 , fullset = " << fullset << endl;
    cout << " data.logDetSubmatrix( fullset ) = " << data.logDetSubmatrix( fullset ) << endl;
#endif
    

    rhs2 += data.logDetSubmatrix( fullset );
#if debug_20140911_1pm_findSegmentationFault
    cout << " tag fw34v4eg5b6 " << endl;
#endif

#if pause_initial_cons
    cout << " & data.logDetSubmatrix( fullset ) = " << data.logDetSubmatrix( fullset ) << " , fullset = " << fullset << " tag 89u2r9823"<< endl;
#endif
    for (int j = 0; j < n; j++){
//        bitset<_CONST_INT_MAX_FACTOR_SIZE> jset(fullset);
//        jset.flip(j);
        bitset<_CONST_INT_MAX_FACTOR_SIZE> jset;
        jset.set(j);
        rhs2 += data.logDetDelta( jset, fullset ) * ( z[j] - 1.0 ) ;
#if 0 // pause_initial_cons
        cout << " & data.logDetDelta( jset, fullset ) = " << data.logDetDelta( jset, fullset ) ;
        cout << ",  jset = " << jset<< endl;
#endif
    }
#if debug_20140911_1pm_findSegmentationFault
    cout << " tag fwf342f " << endl;
#endif

    model.add( f - rhs2 <= 0 );

    
#if debug_20140911_1pm_findSegmentationFault
    cout << "tag f92j0923 \n";
#endif
    
#if pause_initial_cons
    env.out() << "f <= rhs : " << endl;
    env.out() << (f <= rhs) << endl;
    env.out() << "f <= rhs2 : " << endl;
    env.out() << (f <= rhs2) << endl;
    {
        cout << " \n\n\n here is the initial_cons!  input something here: \n\n\n " << endl;
        int temppppl;
        cin >> temppppl;
    }
#endif
    
    if (n < 10)
        env.out() << (f <= rhs) << endl;
    rhs.end();
    if (n < 10)
        env.out() << (f <= rhs2) << endl;
    rhs2.end();
    
    
    
    
};









void SparseRegressionProblem :: build_SR_LowerConstraints (IloModel &model, IloBoolVarArray &z, IloNumVar &g)
{
    
    
    if (!tune.initialConstraintFlag_lower)
        return;

#if title_everywhere
    cout << " SparseRegressionProblem :: build_SR_LowerConstraints tag (345b56ejne5) " << endl;
#endif


#if debugMode_positionTitle
    cout << "  \n * * \n * * build_SR_LowerConstraints. (tag t82c9835) " << endl;
#endif

    // we dont add any constraints by now.
    // implicit constratins are g>= 0
    int n = data.getCols();
    //    double cumulative_violated_amount = 0;
    IloRange cut;
    IloExpr  rhs(env);
    bitset<_CONST_INT_MAX_FACTOR_SIZE> chainsetFront;
    bitset<_CONST_INT_MAX_FACTOR_SIZE> chainsetBack;
    rhs += data.logDetSubmatrix_with_y(chainsetFront);
    for ( int j = 0; j < n; j++) {
        chainsetFront.set( j );
#if debugMode && 0
        {
            cout << " - - - j = " << j << endl;
            cout << "  chainsetFront = " << chainsetFront;
            cout << " , chainsetBack = " << chainsetBack << endl;
            cout << "  data.logDetSubmatrix_with_y(chainsetFront) = " << data.logDetSubmatrix_with_y(chainsetFront);
            cout << " , data.logDetSubmatrix_with_y(chainsetBack) = " << data.logDetSubmatrix_with_y(chainsetBack) << endl;
        }
        {
            int temppp;
            cin >> temppp;
        }
#endif
        rhs += ( data.logDetSubmatrix_with_y(chainsetFront) - data.logDetSubmatrix_with_y(chainsetBack) ) * z[j];
        chainsetBack = chainsetFront;
    };
    model.add( g - rhs >= 0 );
    
    
#if debugMode
    cout << "  \n * * \n * * build_SR_LowerConstraints. tag fq3vrewbwe " << endl;
#endif
    if (n < 10)
        env.out() << (g - rhs >= 0) << endl;

    rhs.end();
    

    
};


void SparseRegressionProblem :: build_SR_Cardinality_Constraints (IloModel &model, IloBoolVarArray &z)
{
//    int k = n/2;
//    IloExpr  rhs(env);
//    for ( int j = 0; j < n; j++) {
//        rhs += z[j];
//    };
////    model.add( rhs <= k);
};








/*
// 
//                  #######  ##     ## ######## ########  ##     ## ########
//                 ##     ## ##     ##    ##    ##     ## ##     ##    ##
//                 ##     ## ##     ##    ##    ##     ## ##     ##    ##
//                 ##     ## ##     ##    ##    ########  ##     ##    ##
//                 ##     ## ##     ##    ##    ##        ##     ##    ##
//                 ##     ## ##     ##    ##    ##        ##     ##    ##
//                  #######   #######     ##    ##         #######     ##
//                 
 */

void SparseRegressionProblem :: output_SR_Solutions( IloCplex &cplex, IloBoolVarArray &z, IloNumVar &f, IloNumVar &g, IloNumVar &mse, TuneParameters &tune )
{
#if title_everywhere
    cout << " SparseRegressionProblem :: output_SR_Solutions tag (fq24bw45b) " << endl;
#endif

#if debugMode_positionTitle
    cout << "  \n * * \n * * output_SR_Solutions. tag f3bil,dr6hnse5 " << endl;
#endif
    
    int n = data.getCols();
    IloNumArray valz(env);
    IloNum valf;
    IloNum valg;
    IloNum valmse;
    cplex.getValues(valz, z);
    valf = cplex.getValue(f);
    valg = cplex.getValue(g);
    valmse = cplex.getValue(mse);
    env.out() << "Solution status = " << cplex.getStatus() << endl;
    env.out() << "Solution value  = " << cplex.getObjValue() << endl;
    env.out() << "# nodes = " << cplex.getNnodes64() << endl;
    env.out() << "val.z  = " << valz << endl;
    env.out() << "valf[u] = " << valf << endl;
    env.out() << "valg[l] = " << valg << endl;
    env.out() << "valmse = " << valmse << endl;
    
    env.out() << "info of cuts: " << endl;
    env.out() << "num_lazy_cut_on_f[u]      =  " << tune.num_lazy_cut_on_f << endl;
    env.out() << "num_lazy_cut_on_g[l]      =  " << tune.num_lazy_cut_on_g << endl;
    env.out() << "num_lazy_cut_on_fandg[u]      =  " << tune.num_lazy_cut_on_fandg << endl;
    env.out() << "num_user_cut_on_f[u]      =  " << tune.num_user_cut_on_f << endl;
    env.out() << "num_user_cut_on_g[l]      =  " << tune.num_user_cut_on_g << endl;
    env.out() << "num_user_cut_on_fandg[u]      =  " << tune.num_user_cut_on_fandg << endl;

    env.out() << "num_lazy_cut           =  " << tune.num_lazy_cut << endl;
    env.out() << "num_user_cut           =  " << tune.num_user_cut << endl;
    env.out() << "num_sbm_lower_cut_on_g[l] =  " << tune.num_sbm_lower_cut_on_g << endl;
    env.out() << "num_sbm_upper_cut_on_f[u] =  " << tune.num_sbm_upper_cut_on_f << endl;
    env.out() << "num_gradient_upper_cut_on_f[u] =  " << tune.num_gradient_upper_cut_on_f << endl;
    env.out() << "num_gradient_upper_cut_on_LMinusU_psd =  " << tune.num_gradient_upper_cut_on_LMinusU_psd << endl;
    env.out() << "num_gradient_upper_cut_on_LMinusU_nsd =  " << tune.num_gradient_upper_cut_on_LMinusU_nsd << endl;
    env.out() << "num_gradient_upper_cut_on_LMinusU_partial =  " << tune.num_gradient_upper_cut_on_LMinusU_partial << endl;
    env.out() << "num_gradient_cut_on_LMinusU_adaptive_psd (!)=  " << tune.num_gradient_cut_on_LMinusU_adaptive_psd << endl;
    env.out() << "num_gradient_cut_on_LMinusU_adaptive_nsd (!)=  " << tune.num_gradient_cut_on_LMinusU_adaptive_nsd << endl;
    
    
//        tune.num_lazy_cut = 0;
//    tune.num_user_cut = 0;
//    
//    tune.num_lazy_cut_on_f = 0;
//    tune.num_lazy_cut_on_g = 0;
//    tune.num_user_cut_on_f = 0;
//    tune.num_user_cut_on_g = 0;
//    tune.num_sbm_lower_cut_on_g = 0;
//    tune.num_sbm_upper_cut_on_f = 0;
//    tune.num_gradient_upper_cut_on_f = 0;

    
    cout << "\n log det:" << std::endl;
    cout << "size = " << data.calculatedValues.size() << std::endl;
    cout << "full.size = " << pow(2,n) << std::endl;
    cout << "bucket_count = " << data.calculatedValues.bucket_count() << std::endl;
    cout << "load_factor = " << data.calculatedValues.load_factor() << std::endl;
    cout << "max_load_factor = " << data.calculatedValues.max_load_factor() << std::endl;
    cout << "\n\n" << std::endl;
    
    
    cout << "\n Node density:" << std::endl;
    cout << "# Node = " << cplex.getNnodes64() << std::endl;
    cout << "# whole space = " << pow(2.0,n) << std::endl;
    cout << "* Node density = " << cplex.getNnodes64() / pow(2.0,n) << std::endl;
    cout << "\n\n" << std::endl;
    
    
    tune.rootGap = (cplex.getObjValue() - tune.rootBestValue) / (1e-10 + abs(cplex.getObjValue()));
    cout << "cplex.getObjValue()  = " << cplex.getObjValue() << ", tune.rootBestValue= " << tune.rootBestValue << endl;
    cout << "rootgap = " << tune.rootGap  << " = " << tune.rootGap*100 << "%" << endl;
    
    cout << "\n"<<cplex.getObjValue()<< "\t" << tune.rootBestValue << "\t" << tune.rootGap*100 << endl;
    cout << "\n\n" << std::endl;
    
    
#if 1
    cout << " scalerLatNSD \n";
    for (int i=0; i<n; i++) {
        cout << data.scalerLatNSD[i] << ", ";
    }
    cout << "\n scalerUatNSD \n";
    for (int i=0; i<n; i++) {
        cout << data.scalerUatNSD[i] << ", ";
    }
    cout << "\n";
#endif
    
    cout << "\n\n" << std::endl;



};























































void SparseRegressionProblem :: generateCuts_total_gradient_adaptive_nsd_copy2( IloEnv &env, IloConstraintArray &cuts2BeAdded, DataModel &data, TuneParameters &tune, IloBoolVarArray &z, IloNumVar &u, IloNumVar &l, IloNumArray valz, IloNum valu, IloNum vall, long nnode, vector<double>  &lastValz)
{
#if test_adap_nsd_cut
    cout << " SparseRegressionProblem :: generateCuts_total_gradient_adaptive_nsd tag (fbw4j7r56hns4e5rtdyc4awe) " << endl;
#endif
    int n = data.getCols();
    
    
    
    double fractional_part = 0;
    for (int i = 0; i < n; i++){
        fractional_part += min( abs(valz[i]), abs(valz[i]-1) );
    }
    if (fractional_part < 1e-4 * n)
        return;
    
    double distFromLastValz = 0;
    for (int i = 0; i < n ; i++){
        ;
        distFromLastValz += pow(lastValz[i] - valz[i], 2);
    }
    cout << " distFromLastValz = " << distFromLastValz << endl;
    if (distFromLastValz < .01*n)
        return;
    
    
    
    double epsilon = .05;
    IloNumArray centralizedValz(env,n);
    for (int i = 0; i < n; i++){
        // find the scale from centralized z.
        // find the gradient from the read z
        if (valz[i] < epsilon )
            centralizedValz[i] = epsilon;
        else if (valz[i] > 1.0-epsilon)
            centralizedValz[i] = 1.0-epsilon;
        else
            centralizedValz[i] = valz[i];
    }
#if test_adap_nsd_cut
    cout << " SparseRegressionProblem :: generateCuts_total_gradient_adaptive_nsd tag (aw4bfyigv4erd) " << endl;
#endif
    
    vector<double> oneMinusXhatOverXhat(n);
    for (int i = 0; i < n; i++)
        oneMinusXhatOverXhat[i] = (1-centralizedValz[i])/centralizedValz[i];
    vector<double> sqrtOneMinusXhatOverXhat(n);
    for (int i = 0; i < n; i++)
        sqrtOneMinusXhatOverXhat[i] = sqrt(oneMinusXhatOverXhat[i]);
    
    vector<vector<double>> XinvLX( data.KUK(sqrtOneMinusXhatOverXhat, data.OriginalLinverse) );
    vector<vector<double>> XinvUX( data.KUK(sqrtOneMinusXhatOverXhat, data.OriginalUinverse) );
    
#if test_adap_nsd_cut
    cout << " SparseRegressionProblem :: generateCuts_total_gradient_adaptive_nsd tag (f34nrgwv4b7na4w) " << endl;
#endif
    
    vector<double> J(n);
    vector<double> K(n);
    for (int i = 0; i < n; i++){
        J[i] = sqrt(1.0/data.OriginalL[i][i]);
        K[i] = sqrt(1.0/data.OriginalU[i][i]);
    }
#if test_adap_nsd_cut || test_adap_nsd_cut_l2
    cout << " tag hrdgerxferbghdrt6nuft ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! " << endl;
    cout << " K = "  << endl;
    for (int i = 0; i < n; i++) cout << K[i] << "\t"; cout << endl;
    cout << " J = "  << endl;
    for (int i = 0; i < n; i++) cout << J[i] << "\t"; cout << endl;
#endif
    
    
    // fit condition 2
    {
        double maxEigOfKUKplusJLJ = data.eig(data.KUKplusJLJ(K,J,data.OriginalU,data.OriginalL))[n-1];
        double scalerForInitialJK = sqrt(2.0/ maxEigOfKUKplusJLJ);
        for (int i = 0; i < n; i++){
            J[i] *= scalerForInitialJK ;
            K[i] *= scalerForInitialJK ;
        };
    }
#if test_adap_nsd_cut || test_adap_nsd_cut_l2
    cout << " tag xtgnhr6dnse4vges ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! , after fit condition 2" << endl;
    cout << " K = "  << endl;
    for (int i = 0; i < n; i++) cout << K[i] << "\t"; cout << endl;
    cout << " J = "  << endl;
    for (int i = 0; i < n; i++) cout << J[i] << "\t"; cout << endl;
#endif
    
    // fit condition 1
    {
        if (data.eig(data.KUKminusJLJ(K, J, data.OriginalU, data.OriginalL))[n-1] > 1e-15){
            double scalerOfK = data.minEigAsqrtinvBAsqrtinv( data.KUK(K,data.OriginalU),data.KUK(J,data.OriginalL) );
            for (int i = 0; i < n; i++)
                K[i] = K[i]*scalerOfK;
        }
    }
#if test_adap_nsd_cut || test_adap_nsd_cut_l2
    cout << " tag fzsefzresvge54bhr6d ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! , after fit condition 21" << endl;
    cout << " K = "  << endl;
    for (int i = 0; i < n; i++) cout << K[i] << "\t"; cout << endl;
    cout << " J = "  << endl;
    for (int i = 0; i < n; i++) cout << J[i] << "\t"; cout << endl;
#endif
    
    // fit condition 2 again
    {
        double maxEigOfKUKplusJLJ = data.eig(data.KUKplusJLJ(K,J,data.OriginalU,data.OriginalL))[n-1];
        double scalerForInitialJK = sqrt(2.0/ maxEigOfKUKplusJLJ);
        for (int i = 0; i < n; i++){
            J[i] *= scalerForInitialJK ;
            K[i] *= scalerForInitialJK ;
        };
    }
    
    
    
#if test_adap_nsd_cut || test_adap_nsd_cut_l2
    cout << " tag 4f35rvgestf3grthrt ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! , after fit condition 212" << endl;
    cout << " K = "  << endl;
    for (int i = 0; i < n; i++) cout << K[i] << "\t"; cout << endl;
    cout << " J = "  << endl;
    for (int i = 0; i < n; i++) cout << J[i] << "\t"; cout << endl;
#endif
    
    
    
    
    
#if test_adap_nsd_cut
    cout << " tag w4ehb5dts4e5t ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! " << endl;
    
    cout << " eig = "  << endl;
    for (int i = 0; i < n; i++) cout << data.eig(data.KUKplusJLJ(K,J,data.OriginalU,data.OriginalL))[i] << "\t"; cout << endl;
    cout << " K = "  << endl;
    for (int i = 0; i < n; i++) cout << K[i] << "\t"; cout << endl;
    cout << " J = "  << endl;
    for (int i = 0; i < n; i++) cout << J[i] << "\t"; cout << endl;
#endif
    
    double deltaTotalObjective = 1e3;
    int maxTotalIterate = 10;
    int iteTotal = 0;
#if test_adap_nsd_cut
    cout << " SparseRegressionProblem :: generateCuts_total_gradient_adaptive_nsd tag (fw3a4v5r65e4awesfx) " << endl;
#endif
    
    
    double totalObjectiveOld;
    {
        vector<double> Ksquare(K.size());
        for (int i = 0; i < n; i++)
            Ksquare[i] = K[i]*K[i];
        vector<double> Jsquare(J.size());
        for (int i = 0; i < n; i++)
            Jsquare[i] = J[i]*J[i];
        totalObjectiveOld = data.logDet_AplusX(Jsquare, XinvLX)-data.logDet_AplusX(Ksquare, XinvUX);
        for(int i =0; i < n; i++){
            totalObjectiveOld += centralizedValz[i]*log(Ksquare[i]) - centralizedValz[i]*log(Jsquare[i]);
        }
    }
#if test_adap_nsd_cut || test_adap_nsd_cut_l2
    cout << " tag f34fawfaw4 ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! " << endl;
    cout << " totalObjectiveOld = " << totalObjectiveOld << endl;
#endif
    
    
    while ( (iteTotal++ < maxTotalIterate) && ( deltaTotalObjective > 1e-1 ) ) {
#if test_adap_nsd_cut
        cout << " tag a3owfj8a4890jfa40 " << endl;
        cout << " iteTotal = " << iteTotal << endl;
        cout << " maxTotalIterate = " << maxTotalIterate << endl;
        cout << " deltaTotalObjective = " << deltaTotalObjective << endl;
#endif
        // move K
        double deltaKObjective = 1e3;
        int iteK = 0;
        int maxKIterate = 10;
        double KstepLength = 5/.618;
#if test_adap_nsd_cut || test_adap_nsd_cut_l2
        cout << " tag aw4v5ge6hrb54dr5  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * " << endl;
        cout << " iteK = " << iteK << endl;
        cout << " maxKIterate = " << maxKIterate << endl;
        cout << " deltaKObjective = " << deltaKObjective << endl;
        cout << " KstepLength = " << KstepLength << endl;
        
#endif
        
        while ( (iteK++ < maxKIterate) && ( deltaKObjective > 1e-1 ) && ( KstepLength > 1e-2 ) ) {
            vector<double> Ksquare(K.size());
            for (int i = 0; i < n; i++)
                Ksquare[i] = K[i]*K[i];
            
            double KobjectiveOld = - data.logDet_AplusX(Ksquare, XinvUX);
            for(int i =0; i < n; i++){
                KobjectiveOld += centralizedValz[i]*log(Ksquare[i]);
            }
            
            
            
            vector<double> gradKsquare( data.gradientLDYplusM( Ksquare, XinvUX));
#if test_adap_nsd_cut
            cout << " tag 4nfjtn4zvcrd " << endl;
            cout << " gradKsquare = "  << endl;
            for (int i = 0; i < n; i++) cout << gradKsquare[i] << "\t"; cout << endl;
#endif
            for (int i = 0; i < n; i++)
                gradKsquare[i] = -gradKsquare[i] + centralizedValz[i]/Ksquare[i] ;
#if test_adap_nsd_cut || test_adap_nsd_cut_l2
            cout << " tag aw34fcgse5b43sz4cf4 " << endl;
            cout << " gradKsquare = "  << endl;
            for (int i = 0; i < n; i++) cout << gradKsquare[i] << "\t"; cout << endl;
#endif
            KstepLength = 5.0/.618;
            bool flagK = true;
            while (flagK){
                KstepLength *= .618;
                if (KstepLength < 1e-2)
                    flagK = false;
                else {
                    vector<double> Knew(n);
                    for (int i=0; i<n; i++)
                        Knew[i] = sqrt(Ksquare[i]+gradKsquare[i]*KstepLength);
                    double minKnew = *min_element(Knew.begin(), Knew.end());
                    
#if test_adap_nsd_cut || test_adap_nsd_cut_l2
                    cout << " tag f4e5gdrbhaw3faw4v " << endl;
                    cout << " KstepLength = " << KstepLength << endl;
                    cout << " minKnew = " << minKnew << endl;
                    cout << " K = "  << endl;
                    for (int i = 0; i < n; i++) cout << K[i] << "\t"; cout << endl;
                    cout << " Knew = "  << endl;
                    for (int i = 0; i < n; i++) cout << Knew[i] << "\t"; cout << endl;
                    cout << " Ksquare = "  << endl;
                    for (int i = 0; i < n; i++) cout << Ksquare[i] << "\t"; cout << endl;
                    cout << " gradKsquare = "  << endl;
                    for (int i = 0; i < n; i++) cout << gradKsquare[i] << "\t"; cout << endl;
#endif
                    
                    
                    if (minKnew < 1e-2)
                        ;
                    else {
                        // scale Knew if needed
                        // vector < vector < double > > KUKminusJLJ( const vector < double > &K , const vector < double > &J , const vector < vector < double > > &U, const vector < vector < double > > &L ) ;
#if test_adap_nsd_cut
                        cout << " tag fw4v4esves5 " << endl;
                        cout << " eig = "  << endl;
                        for (int i = 0; i < n; i++) cout << data.eig(data.KUKminusJLJ(Knew, J, data.OriginalU, data.OriginalL))[i] << "\t"; cout << endl;
#endif
                        if (data.eig(data.KUKminusJLJ(Knew, J, data.OriginalU, data.OriginalL))[n-1] > 1e-15){
                            double scalerOfKnew = data.minEigAsqrtinvBAsqrtinv( data.KUK(Knew,data.OriginalU),data.KUK(J,data.OriginalL) );
                            for (int i = 0; i < n; i++)
                                Knew[i] = Knew[i]*scalerOfKnew;
#if test_adap_nsd_cut || test_adap_nsd_cut_l2
                            cout << " tag fevxd5xd5hx " << endl;
                            cout << " scalerOfKnew = " << scalerOfKnew  << endl;
                            cout << " eig = "  << endl;
                            for (int i = 0; i < n; i++) cout << data.eig(data.KUKminusJLJ(Knew, J, data.OriginalU, data.OriginalL))[i] << "\t"; cout << endl;
                            cout << " Knew = "  << endl;
                            for (int i = 0; i < n; i++) cout << Knew[i] << "\t"; cout << endl;
#endif
                        }
                        //                            // !!!!! to be added: scale Knew
                        //                            ;
                        vector<double> KnewSquare(Knew.size());
                        for (int i = 0; i < n; i++)
                            KnewSquare[i] = Knew[i]*Knew[i];
                        double KobjectiveNew = -data.logDet_AplusX(KnewSquare, XinvUX);
                        for(int i =0; i < n; i++){
                            KobjectiveNew += centralizedValz[i]*log(KnewSquare[i]);
                        }
#if test_adap_nsd_cut || test_adap_nsd_cut_l2
                        
                        cout << " tag 3e5v4cgaew4vae5 " << endl;
                        cout << " KobjectiveNew = " << KobjectiveNew  << endl;
                        cout << " KobjectiveOld = " << KobjectiveOld  << endl;
#endif
                        
                        
                        if ( KobjectiveNew < KobjectiveOld ) {
                            ;
                        } else {
                            flagK = false;
                            deltaKObjective = KobjectiveNew - KobjectiveOld;
#if test_adap_nsd_cut || test_adap_nsd_cut_l2
                            cout << " tag ferbges5nse5ns5en ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! " << endl;
                            cout << " KstepLength = " << KstepLength << endl;
                            cout << " KobjectiveNew = " << KobjectiveNew << endl;
                            cout << " KobjectiveOld = " << KobjectiveOld << endl;
                            cout << " Knew = "  << endl;
                            for (int i = 0; i < n; i++) cout << Knew[i] << "\t"; cout << endl;
                            cout << " K = "  << endl;
                            for (int i = 0; i < n; i++) cout << K[i] << "\t"; cout << endl;
                            cout << " move K ! watch here: %%%%%% %%%%%% %%%%%% %%%%%% %%%%%% %%%%%% %%%%%% %%%%%% %%%%%% %%%%%% %%%%%% %%%%%% %%%%%% %%%%%% %%%%%%  " << endl;
                            
                            {
                                int temp;
                                cin >> temp;
                            }
                            
                            
#endif
                            for (int i=0; i<n; i++)
                                K[i] = Knew[i];
                            KobjectiveOld = KobjectiveNew;
                            
                            
                        }
                    }
                }
            }
        }
#if test_adap_nsd_cut
        cout << " tag faw4ves4ve4, finished choose K ! ! ! ! ! ! " << endl;
        cout << " K = "  << endl;
        for (int i = 0; i < n; i++) cout << K[i] << "\t"; cout << endl;
        cout << " J = "  << endl;
        for (int i = 0; i < n; i++) cout << J[i] << "\t"; cout << endl;
        //        {
        //            int temp;
        //            cin >> temp;
        //        }
#endif
        int maxSIterate = 10;
        double deltaSObjective = 1e3;
        int iteS = 0;
        vector<double> S(n,1.0);
        vector<double> Jinverse(n);
        vector<double> Kinverse(n);
        for (int i = 0; i < n; i++){
            Jinverse[i] = 1.0/J[i];
            Kinverse[i] = 1.0/K[i];
        }
        vector<vector<double>> JinvXinvLXJinv(data.KUK( Jinverse, XinvLX ));
        vector<vector<double>> KinvXinvUXKinv(data.KUK( Kinverse, XinvUX ));
        double SstepLength = 1/.618;
        
        
#if test_adap_nsd_cut
        cout << " tag aw3vaw4vaw4va4w ~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
        cout << " JinvXinvLXJinv = " << endl;
        for (int i = 0; i < n; i++){
            for (int j = 0; j < n; j++)
                cout << JinvXinvLXJinv[i][j] << "\t";
            cout << endl;
        }
        cout << endl;
        
        cout << " KinvXinvUXKinv = " << endl;
        for (int i = 0; i < n; i++){
            for (int j = 0; j < n; j++)
                cout << KinvXinvUXKinv[i][j] << "\t";
            cout << endl;
        }
        cout << endl;
        
        cout << " Jinverse = "  << endl;
        for (int i = 0; i < n; i++) cout << Jinverse[i] << "\t"; cout << endl;
        cout << " Kinverse = "  << endl;
        for (int i = 0; i < n; i++) cout << Kinverse[i] << "\t"; cout << endl;
        
#endif
        
        
        
        
        while ( (iteS++ < maxSIterate) && ( deltaSObjective > 1e-2 ) ) { // && ( SstepLength > 1e-2 ) ) {
#if test_adap_nsd_cut || test_adap_nsd_cut_l2
            cout << " tag fa4v4eab4a " << endl;
            cout << " iteS = " << iteS << endl;
            cout << " maxSIterate = " << maxSIterate << endl;
            cout << " deltaSObjective = " << deltaSObjective << endl;
            cout << " SstepLength = " << SstepLength << endl;
            cout << " S = "  << endl;
            for (int i = 0; i < n; i++) cout << S[i] << "\t"; cout << endl;
#endif
            
            double SobjectiveOld = data.logDet_AplusX(S, JinvXinvLXJinv) - data.logDet_AplusX(S, KinvXinvUXKinv);
            double SobjectiveNextS = SobjectiveOld;
            vector<double> gradS1( data.gradientLDYplusM( S, JinvXinvLXJinv));
            vector<double> gradS2( data.gradientLDYplusM( S, KinvXinvUXKinv));
            vector<double> gradS(n);
            for (int i = 0; i < n; i++)
                gradS[i] = gradS1[i] - gradS2[i] ;
            
            double normGradS = 0;
            for (int i = 0; i < n; i++)
                normGradS += gradS[i] * gradS[i] ;
            normGradS = sqrt(normGradS);
            
            
#if test_adap_nsd_cut
            cout << " tag faw4ves4va4wvtes5bvhr" << endl;
            cout << " S = "  << endl;
            for (int i = 0; i < n; i++) cout << S[i] << "\t"; cout << endl;
            cout << " gradS1 = "  << endl;
            for (int i = 0; i < n; i++) cout << gradS1[i] << "\t"; cout << endl;
            cout << " gradS2 = "  << endl;
            for (int i = 0; i < n; i++) cout << gradS2[i] << "\t"; cout << endl;
            cout << " gradS = "  << endl;
            for (int i = 0; i < n; i++) cout << gradS[i] << "\t"; cout << endl;
            
            //            {
            //                int temp;
            //                cin >> temp;
            //            }
#endif
            
            SstepLength = 100/.618;
            bool flagS = true;
            
            while (flagS){
                SstepLength *= .618;
#if test_adap_nsd_cut
                cout << " normGradS * SstepLength = " << normGradS * SstepLength << endl;
#endif
                if ( ( normGradS * SstepLength < 1e-2*n) || (SstepLength < 1e-2) )
                    flagS = false;
                else {
                    vector<double> Snew(n);
                    for (int i=0; i<n; i++)
                        Snew[i] = S[i]+gradS[i]*SstepLength;
                    double minSnew = *min_element(Snew.begin(), Snew.end());
                    
#if test_adap_nsd_cut || test_adap_nsd_cut_l2
                    cout << " tag fwa4vse4bts35b4  ~ ~ ~ ~ ~ " << endl;
                    cout << " SstepLength = " << SstepLength << endl;
                    cout << " minSnew = " << minSnew << endl;
                    cout << " S = "  << endl;
                    for (int i = 0; i < n; i++) cout << S[i] << "\t"; cout << endl;
                    cout << " Snew = "  << endl;
                    for (int i = 0; i < n; i++) cout << Snew[i] << "\t"; cout << endl;
                    cout << " gradS = "  << endl;
                    for (int i = 0; i < n; i++) cout << gradS[i] << "\t"; cout << endl;
                    
                    
                    cout << " \t \t\t\t\t\t\t\t\t\t\t\t check Snew and S and GradS !!! \n watch here: %%%%%% %%%%%% %%%%% %%% % % % %%%%%% %%%% tag favw4g5be5 %% %%%%%% %%%%%% %%%%%% %%%%%% %%%%%% %%%%%% %%%%%% %%%%%%  " << endl;
                    
                    //                    {
                    //                        int temp;
                    //                        cin >> temp;
                    //                    }
#endif
                    
                    if (minSnew < 1e-2)
                        ;
                    else {
                        vector<double> SnewSqrt(n);
                        for (int i=0; i<n; i++)
                            SnewSqrt[i] = sqrt(Snew[i]);
                        double maxEigSKUKJLJS = data.eig( data.KUK(SnewSqrt, data.KUKplusJLJ(K, J, data.OriginalU, data.OriginalL)) )[n-1]  ;
#if test_adap_nsd_cut
                        cout << " tag fawvf4wbg4eb  ~ ~ ~ ~ ~ " << endl;
                        cout << " SstepLength = " << SstepLength << endl;
                        cout << " maxEigSKUKJLJS = " << maxEigSKUKJLJS << endl;
                        
                        vector<vector<double>> tempM(data.KUKplusJLJ(K, J, data.OriginalU, data.OriginalL));
                        cout << " KUKplusJLJ = " << endl;
                        for (int i = 0; i < n; i++){
                            for (int j = 0; j < n; j++)
                                cout << tempM[i][j] << "\t";
                            cout << endl;
                        }
                        
                        vector<vector<double>> tempM2(data.KUK(SnewSqrt, data.KUKplusJLJ(K, J, data.OriginalU, data.OriginalL)));
                        cout << " S KUKplusJLJ S = " << endl;
                        for (int i = 0; i < n; i++){
                            for (int j = 0; j < n; j++)
                                cout << tempM2[i][j] << "\t";
                            cout << endl;
                        }
                        
                        cout << endl;
                        cout << " eig = "  << endl;
                        for (int i = 0; i < n; i++) cout << data.eig( data.KUK(SnewSqrt, data.KUKplusJLJ(K, J, data.OriginalU, data.OriginalL)) )[i] << "\t"; cout << endl;
                        cout << " S = "  << endl;
                        for (int i = 0; i < n; i++) cout << S[i] << "\t"; cout << endl;
                        cout << " Snew = "  << endl;
                        for (int i = 0; i < n; i++) cout << Snew[i] << "\t"; cout << endl;
                        cout << " SnewSqrt = "  << endl;
                        for (int i = 0; i < n; i++) cout << SnewSqrt[i] << "\t"; cout << endl;
#endif
                        
                        for (int i=0; i<n; i++)
                            Snew[i] = Snew[i]*2/maxEigSKUKJLJS;
                        double SobjectiveNew = data.logDet_AplusX(Snew, JinvXinvLXJinv) - data.logDet_AplusX(Snew, KinvXinvUXKinv);
#if test_adap_nsd_cut|| test_adap_nsd_cut_l2
                        cout << " tag fwefvw4antr6jty   ~ ~ ~ ~ ~ " << endl;
                        
                        cout << " centralizedValz = "  << endl;
                        for (int i = 0; i < n; i++) cout << centralizedValz[i] << "\t"; cout << endl;
                        
                        cout << " SobjectiveNew = " << SobjectiveNew << endl;
                        cout << " SobjectiveOld = " << SobjectiveOld << endl;
                        cout << " Snew = "  << endl;
                        for (int i = 0; i < n; i++) cout << Snew[i] << "\t"; cout << endl;
                        cout << " S = "  << endl;
                        for (int i = 0; i < n; i++) cout << S[i] << "\t"; cout << endl;
                        //                        {
                        //                            int temp;
                        //                            cin >> temp;
                        //                        }
#endif
                        
                        if ( SobjectiveNew < SobjectiveOld ) {
                            ;
                        } else {
                            flagS = false;
#if test_adap_nsd_cut || test_adap_nsd_cut_l2
                            cout << " tag mf7yiftndr5bs ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! " << endl;
                            cout << " SstepLength = " << SstepLength << endl;
                            cout << " SobjectiveOld = " << SobjectiveOld << endl;
                            cout << " SobjectiveNew = " << SobjectiveNew << endl;
                            cout << " S = "  << endl;
                            for (int i = 0; i < n; i++) cout << S[i] << "\t"; cout << endl;
                            cout << " Snew = "  << endl;
                            for (int i = 0; i < n; i++) cout << Snew[i] << "\t"; cout << endl;
                            cout << " gradS = "  << endl;
                            for (int i = 0; i < n; i++) cout << gradS[i] << "\t"; cout << endl;
                            
                            //                            cout << " iteS = " << iteS << endl;
                            //                            cout << " maxSIterate = " << maxSIterate << endl;
                            //                            cout << " deltaSObjective = " << deltaSObjective << endl;
                            //                            cout << " SstepLength = " << SstepLength << endl;
                            cout << " move S !!! watch here: %%%%%% %%%%%% %%%%% %%% % % % %%%%%% %%%% tag ihkmy7fmft %% %%%%%% %%%%%% %%%%%% %%%%%% %%%%%% %%%%%% %%%%%% %%%%%%  " << endl;
                            {
                                int temp;
                                cin >> temp;
                            }
                            
#endif
                            for (int i=0; i<n; i++)
                                S[i] = Snew[i];
                            SobjectiveNextS = SobjectiveNew;
                            deltaSObjective = SobjectiveNextS - SobjectiveOld;
                            SobjectiveOld = SobjectiveNextS;
                        }
                    }
                }
            }
            
            //            deltaSObjective = SobjectiveNew - SobjectiveOld;
            //#if test_adap_nsd_cut || test_adap_nsd_cut_l2
            //            cout << " tag f43vf343b5wnt7 ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! " << endl;
            //            cout << " SobjectiveNew = " << SobjectiveNew << endl;
            //            cout << " SobjectiveOld = " << SobjectiveOld << endl;
            //            cout << " S = "  << endl;
            //            for (int i = 0; i < n; i++) cout << S[i] << "\t"; cout << endl;
            //            cout << " iteS = " << iteS << endl;
            //            cout << " maxSIterate = " << maxSIterate << endl;
            //            cout << " deltaSObjective = " << deltaSObjective << endl;
            //            cout << " SstepLength = " << SstepLength << endl;
            //            {
            //                int temp;
            //                cin >> temp;
            //            }
            //#endif
            //            SobjectiveOld = SobjectiveNew;
        }
        for (int i = 0; i < n; i++){
            double sqrtSi = sqrt(S[i]);
            J[i] = J[i] * sqrtSi;
            K[i] = K[i] * sqrtSi;
        }
        
        {
            vector<double> Ksquare(K.size());
            for (int i = 0; i < n; i++)
                Ksquare[i] = K[i]*K[i];
            vector<double> Jsquare(J.size());
            for (int i = 0; i < n; i++)
                Jsquare[i] = J[i]*J[i];
            double totalObjectiveNew = data.logDet_AplusX(Jsquare, XinvLX)-data.logDet_AplusX(Ksquare, XinvUX);
            for(int i =0; i < n; i++){
                totalObjectiveNew += centralizedValz[i]*log(Ksquare[i]) - centralizedValz[i]*log(Jsquare[i]);
            }
            
            
            deltaTotalObjective = totalObjectiveNew - totalObjectiveOld;
            
#if test_adap_nsd_cut || test_adap_nsd_cut_l2
            cout << " tag fawvfw4cwa, finished choose S ! ! ! ! ! ! " << endl;
            cout << " centralizedValz = "  << endl;
            for (int i = 0; i < n; i++) cout << centralizedValz[i] << "\t"; cout << endl;
            cout << " real z = "  << endl;
            for (int i = 0; i < n; i++) cout << valz[i] << "\t"; cout << endl;
            cout << " S = "  << endl;
            for (int i = 0; i < n; i++) cout << S[i] << "\t"; cout << endl;
            cout << " K = "  << endl;
            for (int i = 0; i < n; i++) cout << K[i] << "\t"; cout << endl;
            cout << " J = "  << endl;
            for (int i = 0; i < n; i++) cout << J[i] << "\t"; cout << endl;
            
            
            cout << " totalObjectiveNew = " << totalObjectiveNew  << endl;
            cout << " totalObjectiveOld = " << totalObjectiveOld  << endl;
            cout << " deltaTotalObjective = " << deltaTotalObjective  << endl;
            cout << " iteTotal = " << iteTotal  << endl;
            
            
            {
                int temp;
                cin >> temp;
            }
#endif
            totalObjectiveOld = totalObjectiveNew;
#if test_adap_nsd_cut || test_adap_nsd_cut_l2
            cout << " tag fawefwafw443 ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! " << endl;
            cout << " totalObjectiveOld = " << totalObjectiveOld << ", iteTotal = " << iteTotal << endl;
#endif
            
#if test_adap_nsd_cut || test_adap_nsd_cut_l2
            cout << " tag zsevfaw4bdr6nut7  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * " << endl;
            cout << " iteTotal = " << iteTotal << endl;
            cout << " maxTotalIterate = " << maxTotalIterate << endl;
            cout << " deltaTotalObjective = " << deltaTotalObjective << endl;
            cout << " condition:     while ( (iteTotal++ < maxTotalIterate) && ( deltaTotalObjective > 1e-2 ) )  " << endl;
#endif
            
            
            
        }
        
        
        
        
    }
    
    
    
    
    
    
    
#if test_adap_nsd_cut || test_adap_nsd_cut_l2
    cout << " tag f4wbtes5br65n tujtf ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !  ready to add adaptive nsd cut" << endl;
    cout << " centralizedValz = "  << endl;
    for (int i = 0; i < n; i++) cout << centralizedValz[i] << "\t"; cout << endl;
    cout << " real z = "  << endl;
    for (int i = 0; i < n; i++) cout << valz[i] << "\t"; cout << endl;
    cout << " K = "  << endl;
    for (int i = 0; i < n; i++) cout << K[i] << "\t"; cout << endl;
    cout << " J = "  << endl;
    for (int i = 0; i < n; i++) cout << J[i] << "\t"; cout << endl;
    {
        int temp;
        cin >> temp;
    }
#endif
    
    
    // from here, we can add cut according to J and K
    
    
    IloRange cut;
    IloExpr rhs(env);
    try {
        vector < vector < double > > scaledL = data.scaleBothSide(data.OriginalL, J);;
        vector < vector < double > > scaledU = data.scaleBothSide(data.OriginalU, K);;
        
#if test_adap_nsd_cut && 0
        cout << " \n\n\n\n tag (7jue5bh63w4v4wsger) " << endl ;
        cout << " data.OriginalL = " << endl ;
        for (int i = 0; i < n; i++){
            for (int j = 0; j < n; j++)
                cout << data.OriginalL[i][j] << "\t" ;
            cout << endl;
        }
        cout << " data.OriginalU = " << endl ;
        for (int i = 0; i < n; i++){
            for (int j = 0; j < n; j++)
                cout << data.OriginalU[i][j] << "\t" ;
            cout << endl;
        }
        cout << " scaledL = " << endl ;
        for (int i = 0; i < n; i++){
            for (int j = 0; j < n; j++)
                cout << scaledL[i][j] << "\t" ;
            cout << endl;
        }
        cout << " scaledU = " << endl ;
        for (int i = 0; i < n; i++){
            for (int j = 0; j < n; j++)
                cout << scaledU[i][j] << "\t" ;
            cout << endl;
        }
        cout << " valz =  " ;
        for (int i = 0; i < n; i++)
            cout << valz[i] << "\t" ;
        cout << endl;
        cout << " J =  " ;
        for (int i = 0; i < n; i++)
            cout << J[i] << "\t" ;
        cout << endl;
        cout << " K =  " ;
        for (int i = 0; i < n; i++)
            cout << K[i] << "\t" ;
        cout << endl;
        
        {
            int temp;
            cin >> temp;
        }
#endif
        
        
        
        double lpart = data.logDet_IplusDiagzW( valz, scaledL )  + log(data.Sigma[n][n]) ;
        double upart = data.logDet_IplusDiagzW( valz, scaledU ) ;
        for (int i = 0; i < n; i++){
            lpart -= 2*log(J[i])*valz[i];
            upart -= 2*log(K[i])*valz[i];
        }
        
        double c = lpart - upart;
        
        /* test if current solution violates the cutting planes. if yes, add the cut */
        double criterior2 = (vall - valu - c);
#if test_adap_nsd_cut || test_adap_nsd_cut_l2
        cout << "   * adaptive nsd cut * * adaptive nsd cut * * adaptive nsd cut * * adaptive nsd cut * * adaptive nsd cut * * adaptive nsd cut * * adaptive nsd cut * * adaptive nsd cut * * adaptive nsd cut *  " << endl;
        cout << "c = " << c << endl;
        cout << "lpart = " << lpart << ", upart = " << upart << ", lpart-upart = " << lpart-upart << endl;
        cout << "vall = " << vall << ", valu = " << valu << ", vall-valu = " << vall-valu << endl;
        cout << " ** criterior2 = (vall - valu - c) = " << criterior2 << endl;
        cout << " valz = "  << valz << endl;
#endif
        
        // rhs should be 0 now;
        if ( ( (nnode < 1000) && ( criterior2 < -1e-2 ) )  || ( (nnode > 1000) && ( criterior2 < -1e-1 ) ) ) {
            rhs += lpart - upart;
            vector<double> grad( data.uMinusLGradient( valz, scaledL, scaledU ) );
            for (int i = 0; i < n; i++){
                grad[i] -= 2*log(J[i]);
                grad[i] += 2*log(K[i]);
            }
            
            for (int i = 0; i < n; i++){
                rhs += grad[i] * ( z[i] - valz[i] );
            }
            
            cut = ( l - u - rhs >= 0 );
#if test_adap_nsd_cut || test_adap_nsd_cut_l2
            cout << " cut = "  << cut << endl;
            {
                int temp;
                cin >> temp;
            }
#endif
            
            cuts2BeAdded.add(cut);
            tune.num_gradient_cut_on_LMinusU_adaptive_nsd++;
        }
        rhs.end();
    } catch (...) {
        rhs.end();
        throw;
    };
    
    
    
    
    
};


void SparseRegressionProblem :: generateCuts_total_gradient_adaptive_psd_copy2( IloEnv &env, IloConstraintArray &cuts2BeAdded, DataModel &data, TuneParameters &tune, IloBoolVarArray &z, IloNumVar &u, IloNumVar &l, IloNumArray valz, IloNum valu, IloNum vall, long nnode, vector<double>  &lastValz)
{
#if test_adap_psd_cut
    cout << " SparseRegressionProblem :: generateCuts_total_gradient_adaptive_psd tag (gw3ga4egw345q) " << endl;
#endif
    int n = data.getCols();
    
    
    
    double fractional_part = 0;
    for (int i = 0; i < n; i++){
        fractional_part += min( abs(valz[i]), abs(valz[i]-1) );
    }
    if (fractional_part < 1e-4 * n)
        return;
    
    double distFromLastValz = 0;
    for (int i = 0; i < n ; i++){
        ;
        distFromLastValz += pow(lastValz[i] - valz[i], 2);
        
    }
    cout << " distFromLastValz = " << distFromLastValz << endl;
    if (distFromLastValz < .01*n)
        return;
    
    
    
    
    
    
    
    double epsilon = .1;
    IloNumArray centralizedValz(env,n);
    for (int i = 0; i < n; i++){
        // find the scale from centralized z.
        // find the gradient from the read z
        if (valz[i] < epsilon )
            centralizedValz[i] = epsilon;
        else if (valz[i] > 1.0-epsilon)
            centralizedValz[i] = 1.0-epsilon;
        else
            centralizedValz[i] = valz[i];
    }
#if test_adap_psd_cut
    cout << " SparseRegressionProblem :: generateCuts_total_gradient_nsd tag (afwvges54b) " << endl;
#endif
    
    
    vector<double> oneMinusXhatOverXhat(n);
    for (int i = 0; i < n; i++)
        oneMinusXhatOverXhat[i] = (1-centralizedValz[i])/centralizedValz[i];
    vector<double> sqrtOneMinusXhatOverXhat(n);
    for (int i = 0; i < n; i++)
        sqrtOneMinusXhatOverXhat[i] = sqrt(oneMinusXhatOverXhat[i]);
    vector<double> y(n);
    
    int ITSTEP = 3;  // must be >= 1
    double STEPLENGTH = 1;
    // obj = LD(Y + L) - LD(Y + U)
#if test_adap_psd_cut
    cout << " SparseRegressionProblem :: generateCuts_total_gradient_nsd tag (aw4bgrd65bny68ky) " << endl;
#endif
    for (int ite = 0; ite <= ITSTEP; ite++) {
        // set up
        if (ite == 0) {
            for (int i = 0; i < n; i++)
                y[i] = oneMinusXhatOverXhat[i];
        }
        vector<double> sqrty(n);
        for (int i = 0; i < n; i++)
            sqrty[i] = sqrt(y[i]);
        // shrink y;
        
        
        vector< vector< double > > LHS = data.OriginalInvOfInvSum;
#if test_adap_psd_cut
        cout << " tag (faw4fv5e) " << endl;
        cout << " sqrtOneMinusXhatOverXhat =  " ;
        for (int i = 0; i < n; i++)
            cout << sqrtOneMinusXhatOverXhat[i] << "\t" ;
        cout << endl;
        cout << " sqrty = " ;
        for (int i = 0; i < n; i++)
            cout << sqrty[i] << "\t" ;
        cout << endl;
        cout << " LHS = " ;
        for (int i = 0; i < n; i++){
            for (int j = 0; j < n; j++)
                cout << LHS[i][j] << "\t" ;
            cout << endl;
        }
        cout << endl;
#endif
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                LHS[i][j] *= sqrtOneMinusXhatOverXhat[i] * sqrtOneMinusXhatOverXhat[j] / sqrty[i] / sqrty[j];
        vector<double> e = data.eig(LHS) ;
        
#if test_adap_psd_cut
        cout << " tag (fq34q34g34) " << endl;
        cout << " eigen of LHS =  " << endl;
        for (int i = 0; i < n; i++)
            cout << e[i] << "\t" ;
        cout << endl;
#endif
        for (int i = 0; i < n; i++)
            y[i] /= e[n-1];
        
        // move y
        if (ite == ITSTEP)
            break;
        vector<double> grad = data.gradientLDYplusM(y, data.OriginalL );
        vector<double> gradU = data.gradientLDYplusM(y, data.OriginalU );
        for (int i = 0; i < n; i++)
            grad[i] -= gradU[i];
#if test_adap_psd_cut
        cout << " tag (gq3ga4wga4) " << endl;
        cout << " y =  " ;
        for (int i = 0; i < n; i++)
            cout << y[i] << "\t" ;
        cout << endl;
        cout << " grad =  " ;
        for (int i = 0; i < n; i++)
            cout << grad[i] << "\t" ;
        cout << endl;
#endif
        for (int i = 0; i < n; i++)
            y[i] += STEPLENGTH * grad[i];
    }
    
    vector<double> scale(n);
    for (int i = 0; i < n; i++) {
        scale[i] = sqrt( oneMinusXhatOverXhat[i] / y[i] );
    }
    
    
#if test_adap_psd_cut
    cout << " tag (faw4q3b5sen) " << endl;
    cout << " y =  " ;
    for (int i = 0; i < n; i++)
        cout << y[i] << "\t" ;
    cout << endl;
    cout << " scale =  " ;
    for (int i = 0; i < n; i++)
        cout << scale[i] << "\t" ;
    cout << endl;
#endif
    
    
    /*
     start to add cut !!!
     */
    IloRange cut;
    IloExpr rhs(env);
    try {
        vector < vector < double > > scaledL = data.scaleBothSide(data.OriginalL, scale);;
        vector < vector < double > > scaledU = data.scaleBothSide(data.OriginalU, scale);;
        
        
        
#if test_adap_psd_cut
        cout << " \n\n\n\n tag (faw4zfg567tn65he4s) " << endl ;
        cout << " data.OriginalL = " << endl ;
        for (int i = 0; i < n; i++){
            for (int j = 0; j < n; j++)
                cout << data.OriginalL[i][j] << "\t" ;
            cout << endl;
        }
        cout << " data.OriginalU = " << endl ;
        for (int i = 0; i < n; i++){
            for (int j = 0; j < n; j++)
                cout << data.OriginalU[i][j] << "\t" ;
            cout << endl;
        }
        cout << " scaledL = " << endl ;
        for (int i = 0; i < n; i++){
            for (int j = 0; j < n; j++)
                cout << scaledL[i][j] << "\t" ;
            cout << endl;
        }
        cout << " scaledU = " << endl ;
        for (int i = 0; i < n; i++){
            for (int j = 0; j < n; j++)
                cout << scaledU[i][j] << "\t" ;
            cout << endl;
        }
        cout << " valz =  " ;
        for (int i = 0; i < n; i++)
            cout << valz[i] << "\t" ;
        cout << endl;
        cout << " scale =  " ;
        for (int i = 0; i < n; i++)
            cout << scale[i] << "\t" ;
        cout << endl;
        
        
        
        {
            int temp;
            cin >> temp;
        }
#endif
        
        
        
        double lpart = data.logDet_IplusDiagzW( valz, scaledL )  + log(data.Sigma[n][n]) ;
        double upart = data.logDet_IplusDiagzW( valz, scaledU ) ;
        double c = lpart - upart;
        
        /* test if current solution violates the cutting planes. if yes, add the cut */
        double criterior2 = (vall - valu - c);
#if test_adap_psd_cut
        cout << "   * adaptive psd cut * * adaptive psd cut * * adaptive psd cut * * adaptive psd cut * * adaptive psd cut * * adaptive psd cut * * adaptive psd cut * * adaptive psd cut * * adaptive psd cut *  " << endl;
        cout << "c = " << c << endl;
        cout << "lpart = " << lpart << ", upart = " << upart << ", lpart-upart = " << lpart-upart << endl;
        cout << "vall = " << vall << ", valu = " << valu << ", vall-valu = " << vall-valu << endl;
        cout << " ** criterior2 = (vall - valu - c) = " << criterior2 << endl;
        cout << " valz = "  << valz << endl;
#endif
        
        // rhs should be 0 now;
        if ( ( (nnode < 1000) && ( criterior2 < -1e-2 ) )  || ( (nnode > 1000) && ( criterior2 < -1e-1 ) ) ) {
            rhs += lpart - upart;
            vector<double> grad( data.uMinusLGradient( valz, scaledL, scaledU ) );
            for (int i = 0; i < n; i++){
                rhs += grad[i] * ( z[i] - valz[i] );
            }
            
            cut = ( l - u - rhs >= 0 );
#if test_adap_psd_cut
            cout << " cut = "  << cut << endl;
#endif
            
            cuts2BeAdded.add(cut);
            tune.num_gradient_cut_on_LMinusU_adaptive_psd++;
        }
        
        
        rhs.end();
        
    } catch (...) {
        rhs.end();
        
        throw;
    };
    
};





void SparseRegressionProblem :: generateCuts_total_gradient_adaptive_psd_copy20141031_0138pm( IloEnv &env, IloConstraintArray &cuts2BeAdded, DataModel &data, TuneParameters &tune, IloBoolVarArray &z, IloNumVar &u, IloNumVar &l, IloNumArray valz, IloNum valu, IloNum vall, long nnode, vector<double>  &lastValz)
{
#if test_adap_psd_cut
    cout << " SparseRegressionProblem :: generateCuts_total_gradient_adaptive_psd tag (gw3ga4egw345q) " << endl;
#endif
    int n = data.getCols();
    
    
    
    double fractional_part = 0;
    for (int i = 0; i < n; i++){
        fractional_part += min( abs(valz[i]), abs(valz[i]-1) );
    }
    if (fractional_part < 1e-4 * n)
        return;
    
    double distFromLastValz = 0;
    for (int i = 0; i < n ; i++){
        ;
        distFromLastValz += pow(lastValz[i] - valz[i], 2);
        
    }
    cout << " distFromLastValz = " << distFromLastValz << endl;
    if (distFromLastValz < .01*n)
        return;
    
    
    
    
    double epsilon = .05;
    IloNumArray centralizedValz(env,n);
    for (int i = 0; i < n; i++){
        // find the scale from centralized z.
        // find the gradient from the read z
        if (valz[i] < epsilon )
            centralizedValz[i] = epsilon;
        else if (valz[i] > 1.0-epsilon)
            centralizedValz[i] = 1.0-epsilon;
        else
            centralizedValz[i] = valz[i];
    }
#if test_adap_psd_cut
    cout << " SparseRegressionProblem :: generateCuts_total_gradient_psd tag (ykftjft6udtjth) " << endl;
#endif
    
    
    vector<double> oneMinusXhatOverXhat(n);
    for (int i = 0; i < n; i++)
        oneMinusXhatOverXhat[i] = (1-centralizedValz[i])/centralizedValz[i];
    vector<double> sqrtOneMinusXhatOverXhat(n);
    for (int i = 0; i < n; i++)
        sqrtOneMinusXhatOverXhat[i] = sqrt(oneMinusXhatOverXhat[i]);
    vector< vector< double > > ConstraintLHS = data.OriginalInvOfInvSum;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            ConstraintLHS[i][j] *= sqrtOneMinusXhatOverXhat[i] * sqrtOneMinusXhatOverXhat[j];
    
    
    
    
    
    vector<double> y(n);
    
    for (int i = 0; i < n; i++)
        y[i] = oneMinusXhatOverXhat[i];
    vector<double> sqrty(n);
    for (int i = 0; i < n; i++)
        sqrty[i] = sqrt(y[i]);
    
    
    
    
    vector< vector< double > > C = data.OriginalInvOfInvSum;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            // 要不要乘一个2 ? 不用 he和
            C[i][j] *= sqrtOneMinusXhatOverXhat[i] * sqrtOneMinusXhatOverXhat[j]; // C = X% * 2 inv[ invU + invL] * X%
    
    vector< vector< double > > LHS = C;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            LHS[i][j] = C[i][j] / sqrty[i] / sqrty[j];
    vector<double> e = data.eig(LHS) ;
    for (int i = 0; i < n; i++)
        y[i] *= e[0]/2.0;
    // finished y initialization
    
    
#if test_adap_psd_cut_level2 && 0
    cout << " tag zsefewfwagrwe ";
    cout << " y =  " ;
    for (int i = 0; i < n; i++)
        cout << y[i] << "\t" ;
    cout << endl;
    cout << " just initialize y, current obj =  " << data.logDet_AplusX( y , data.OriginalL ) - data.logDet_AplusX( y , data.OriginalU ) + log(data.Sigma[n][n])<< endl;
    { int temp; cin >> temp; }
#endif
    
#if test_adap_psd_cut_level3
    double TEST_pure_totalObjectiveNew;
    double TEST_pure_totalObjectiveOld;
#endif
    
    
    double deltaTotalObjective = 1e3;
    //    int maxTotalIterate = 5;
    //    int iteTotal = 0;
#if test_adap_nsd_cut
    cout << " SparseRegressionProblem :: generateCuts_total_gradient_adaptive_psd tag (y8tnrd6ybse5ta4) " << endl;
#endif
    double mu = 1;
    
    vector<double> negativey(n);
    for (int i = 0; i < n; i++)
        negativey[i] = -(y[i]);
    double totalObjectiveOld = data.logDet_AplusX( y , data.OriginalL ) - data.logDet_AplusX( y , data.OriginalU ) + log(data.Sigma[n][n]);
#if test_adap_psd_cut_level3
    TEST_pure_totalObjectiveOld =totalObjectiveOld;
#endif
    totalObjectiveOld += mu * data.logDet_AplusX( negativey, C );
    for (int i = 0; i < n; i++)
        totalObjectiveOld += mu * log(y[i]);
    
    
    for (; mu > 1e-6; mu/=10.0) {
#if test_adap_psd_cut_level3
        cout << " SparseRegressionProblem :: generateCuts_total_gradient_adaptive_psd tag (rtsynse5ymt7i6.) " << endl;
        cout << " mu = " << mu << endl;
        {
            int temp;
            cin >> temp;
        }
#endif
        bool flagFindOptimalinMu = false;
        while (!flagFindOptimalinMu) {
            vector<double> direction = data.newtonDirectionPSDSeparationProblem(y, data.OriginalL, data.OriginalU, C, mu);
#if test_adap_psd_cut_level3
            cout << " SparseRegressionProblem :: generateCuts_total_gradient_adaptive_psd tag (fwbr5xnydr6m.) " << endl;
            cout << " y = [ " ;
            for (int i = 0; i < n; i++)
                cout << y[i] << ", " ;
            cout << " ] " << endl;
            cout << " direction = [ ";
            for (int i = 0; i < n; i++)
                cout << direction[i] << ", " ;
            cout << " ] " << endl;
            
#endif
            double normDirection = 0;
            for (int i = 0; i < n; i++)
                normDirection += direction[i]*direction[i];
            normDirection = sqrt(normDirection);
#if test_adap_psd_cut_level3
            cout << " SparseRegressionProblem :: generateCuts_total_gradient_adaptive_psd tag (fawf aerg nhrtm.) " << endl;
            cout << " mu = " << mu << endl;
            cout << " normDirection = " <<normDirection <<  endl;
            cout << " 1e-1*mu*n = " << 1e-1*mu*n << endl;
            cout << " max(1e-1*mu, 1e-3)*n = " << max(1e-1*mu, 1e-3)*n << endl;
            cout << " direction = [ ";
            for (int i = 0; i < n; i++)
                cout << direction[i] << ", " ;
            cout << " ] " << endl;
            {
                int temp;
                cin >> temp;
            }
#endif
            if (normDirection < max(1e-1*mu, 1e-3)*n) {
                flagFindOptimalinMu = true;
            } else {
                double YstepLength = 1;
                double YBestStepLength = 1;
                vector<double> newy(n);
                bool flag = true;
                while ( flag && (YstepLength < 1000) ){
                    YstepLength *= 3;
#if test_adap_psd_cut_level3
                    cout << " SparseRegressionProblem :: generateCuts_total_gradient_adaptive_psd tag (fawef4be5jytu) " << endl;
                    cout << " YstepLength = " << YstepLength << endl;
#endif
                    for (int i = 0; i < n; i++){
                        newy[i] = y[i] + YstepLength * direction[i];
                        if ( newy[i] < 1e-5 ) {
                            flag = false;
                            break;
                        }
                    }
#if test_adap_psd_cut_level3
                    cout << " newy = [ " ;
                    for (int i = 0; i < n; i++)
                        cout << newy[i] << ", " ;
                    cout << " ] " << endl;
                    
                    cout << " flag = " << flag  << " , jfwajf0awjf0 "<< endl;
#endif
                    if (flag) {
#if test_adap_psd_cut_level3
                        cout << " SparseRegressionProblem :: generateCuts_total_gradient_adaptive_psd tag (fawfa4ege5g3q) " << endl;
                        
#endif
                        vector<vector<double>> matrix(n, vector<double>(n));
                        
#if test_adap_psd_cut_level3
                        cout << " SparseRegressionProblem :: generateCuts_total_gradient_adaptive_psd tag (awf4ewge5he5) " << endl;
                        
                        cout << " matrix = " << endl;
                        for (int i = 0; i < n; i++){
                            for (int j = 0; j < n; j++)
                                cout <<matrix[i][j] << ", " ;
                            cout << endl;
                        }
                        
#endif
                        for (int i = 0; i < n; i++)
                            for (int j = 0; j < n; j++)
                                if (i == j)
                                    matrix[i][i] = C[i][i]-newy[i];
                                else
                                    matrix[i][j] = C[i][j];
#if test_adap_psd_cut_level3
                        cout << " SparseRegressionProblem :: generateCuts_total_gradient_adaptive_psd tag (wgaew4geagesh) " << endl;
                        
                        cout << " matrix = " << endl;
                        for (int i = 0; i < n; i++){
                            for (int j = 0; j < n; j++)
                                cout <<matrix[i][j] << ", " ;
                            cout << endl;
                        }
                        
                        cout << " C = " << endl;
                        for (int i = 0; i < n; i++){
                            for (int j = 0; j < n; j++)
                                cout <<C[i][j] << ", " ;
                            cout << endl;
                        }
                        {
                            cout << " ! ! ! ! ~ ~ ~ ~ ! ! ! ! ~ ~ ~ ~ ! ! ! ! ~ ~ ~ ~ ! ! ! ! ~ ~ ~ ~ ! ! ! ! ~ ~ ~ ~ ! ! ! ! ~ ~ ~ ~ " << endl;
                            int temp;
                            cin >> temp;
                        }
                        
                        cout << " newy = [ ";
                        for (int i = 0; i < n; i++)
                            cout << newy[i] << ", " ;
                        cout << " ] " << endl;
                        
                        
                        cout << " data.eig(matrix)[0] = [ ";
                        for (int i = 0; i < n; i++)
                            cout << data.eig(matrix)[i] << ", " ;
                        cout << " ] " << endl;
#endif
                        
                        if ( data.eig(matrix)[0] < 1e-12 )
                            flag = false; // otherwise this is a feasible solution
                        
                    }
#if test_adap_psd_cut_level3
                    cout << " flag = " << flag ;
#endif
                    if (flag) {
                        vector<double> negativey(n);
                        for (int i = 0; i < n; i++)
                            negativey[i] = - newy[i] ;
                        double totalObjectiveNew;
                        totalObjectiveNew = data.logDet_AplusX( newy , data.OriginalL ) - data.logDet_AplusX( newy , data.OriginalU ) + log(data.Sigma[n][n]);
#if test_adap_psd_cut_level3
                        TEST_pure_totalObjectiveNew = totalObjectiveNew;
#endif
                        totalObjectiveNew += mu * data.logDet_AplusX( negativey, C );
                        for (int i = 0; i < n; i++)
                            totalObjectiveNew += mu * log(newy[i]);
                        
#if test_adap_psd_cut_level3
                        cout << " SparseRegressionProblem :: generateCuts_total_gradient_adaptive_psd tag (fq4wfvw4stjntuyk) " << endl;
                        cout << " totalObjectiveNew = " << totalObjectiveNew;
                        cout << ", totalObjectiveOld = " << totalObjectiveOld <<  endl;
                        
#if test_adap_psd_cut_level3
                        cout << " TEST_pure_totalObjectiveNew = " << TEST_pure_totalObjectiveNew;
                        cout << ", TEST_pure_totalObjectiveOld = " << TEST_pure_totalObjectiveOld <<  endl;
#endif
                        
                        
                        
                        cout << " YstepLength = " << YstepLength << ", mu = " << mu << endl;
                        cout << " newy = [ " ;
                        for (int i = 0; i < n; i++)
                            cout << newy[i] << ", " ;
                        cout << " ] " << endl;
#endif
                        
                        if ( totalObjectiveNew - totalObjectiveOld < 0.01*YstepLength)
                            flag = false;
                        else
                            YBestStepLength = YstepLength;
                    }
                };
                flag = true;
                YstepLength = YBestStepLength / .5;
                while (flag && (YstepLength > 1e-3 ) ){
                    bool subFlag = true;
                    YstepLength *= .5;
                    for (int i = 0; i < n; i++){
                        newy[i] = y[i] + YstepLength * direction[i];
                        if ( newy[i] < 1e-5 ) {
                            subFlag = false;
                            break;
                        }
                    }
                    if (subFlag) {
                        vector<vector<double>> matrix(n, vector<double>(n));
                        for (int i = 0; i < n; i++)
                            for (int j = 0; j < n; j++)
                                if (i == j)
                                    matrix[i][i] = C[i][i]-newy[i];
                                else
                                    matrix[i][j] = C[i][j];
                        if ( data.eig(matrix)[0] < 1e-12 )
                            subFlag = false; // otherwise this is a feasible solution
                        
                    }
                    if (subFlag) {
                        vector<double> negativey(n);
                        for (int i = 0; i < n; i++)
                            negativey[i] = - newy[i] ;
                        double totalObjectiveNew;
                        totalObjectiveNew = data.logDet_AplusX( newy , data.OriginalL ) - data.logDet_AplusX( newy , data.OriginalU ) + log(data.Sigma[n][n]);
#if test_adap_psd_cut_level3
                        TEST_pure_totalObjectiveNew = totalObjectiveNew;
#endif
                        totalObjectiveNew += mu * data.logDet_AplusX( negativey, C );
                        for (int i = 0; i < n; i++)
                            totalObjectiveNew += mu * log(newy[i]);
                        
#if test_adap_psd_cut_level3
                        cout << " SparseRegressionProblem :: generateCuts_total_gradient_adaptive_psd tag (wafwa4g) " << endl;
                        cout << " totalObjectiveNew = " << totalObjectiveNew;
                        cout << ", totalObjectiveOld = " << totalObjectiveOld <<  endl;
                        cout << " TEST_pure_totalObjectiveNew = " << TEST_pure_totalObjectiveNew;
                        cout << ", TEST_pure_totalObjectiveOld = " << TEST_pure_totalObjectiveOld <<  endl;
                        
                        cout << " YstepLength = " << YstepLength << endl;
#endif
                        
                        if ( totalObjectiveNew - totalObjectiveOld < 0.01 * YstepLength){
                            ;
                        } else {
                            YBestStepLength = YstepLength;
                            
                            totalObjectiveOld = totalObjectiveNew;
#if test_adap_psd_cut_level3
                            TEST_pure_totalObjectiveOld = totalObjectiveOld;
#endif
                            
                            flag = false;
                        }
                    }
                };
                
                if (flag) {
                    // not need to update y
                    flagFindOptimalinMu = true;
                    
                } else {
#if test_adap_psd_cut_level3
                    cout << " SparseRegressionProblem :: generateCuts_total_gradient_adaptive_psd tag ( fawgw35hrwhr ) " << endl;
                    cout << " update ! ! ! ! ! ! ! ! !! ! ! ! ! ! ! ! ! ! ! y[i]  ! ! ! ! ! ! ! ! !! ! ! ! ! ! ! ! ! ! !  ! ! ! ! ! ! ! ! !! ! ! ! ! ! ! ! ! ! !  ! ! ! ! ! ! ! ! !! ! ! ! ! ! ! ! ! ! !  "  << endl;
                    cout << "YstepLength = " << YstepLength << endl;
                    cout << " y = [ " ;
                    for (int i = 0; i < n; i++)
                        cout << y[i] << ", " ;
                    cout << " ] " << endl;
                    
#endif
                    for (int i = 0; i < n; i++)
                        y[i] = y[i] + YstepLength * direction[i];
                    
#if test_adap_psd_cut_level3
                    
                    cout << " y = [ " ;
                    for (int i = 0; i < n; i++)
                        cout << y[i] << ", " ;
                    cout << " ] " << endl;
#endif
                    
                }
            }
        }
    }
    
    //    {
    //        int ITSTEP = 3;  // must be >= 1
    //        double STEPLENGTH = 1;
    //        // obj = LD(Y + L) - LD(Y + U)
    //#if test_adap_psd_cut
    //        cout << " SparseRegressionProblem :: generateCuts_total_gradient_psd tag (aw4bgrd65bny68ky) " << endl;
    //#endif
    //        for (int ite = 0; ite <= ITSTEP; ite++) {
    //            // set up
    //            if (ite == 0) {
    //                for (int i = 0; i < n; i++)
    //                    y[i] = oneMinusXhatOverXhat[i];
    //
    //            }
    //            vector<double> sqrty(n);
    //            for (int i = 0; i < n; i++)
    //                sqrty[i] = sqrt(y[i]);
    //            // shrink y;
    //
    //
    //            vector< vector< double > > LHS = data.OriginalInvOfInvSum;
    //            for (int i = 0; i < n; i++)
    //                for (int j = 0; j < n; j++)
    //                    // 要不要乘一个2 ? 不用 he和
    //                    LHS[i][j] *= sqrtOneMinusXhatOverXhat[i] * sqrtOneMinusXhatOverXhat[j] / sqrty[i] / sqrty[j];
    //            vector<double> e = data.eig(LHS) ;
    //#if test_adap_psd_cut_level2
    //            cout << " tag fawgeragerhetrher ";
    //            cout << " e =  " ;
    //            for (int i = 0; i < n; i++)
    //                cout << e[i] << "\t" ;
    //            cout << endl;
    //            { int temp; cin >> temp; }
    //#endif
    //
    //
    //#if test_adap_psd_cut
    //            cout << " tag (fq34q34g34) " << endl;
    //            cout << " eigen of LHS =  " << endl;
    //            for (int i = 0; i < n; i++)
    //                cout << e[i] << "\t" ;
    //            cout << endl;
    //#endif
    //            for (int i = 0; i < n; i++)
    //                //            y[i] /= e[n-1];
    //                y[i] *= e[0];
    //#if test_adap_psd_cut_level2
    //            cout << " tag fawgeragerhetrher ";
    //            cout << " e =  " ;
    //            for (int i = 0; i < n; i++)
    //                cout << e[i] << "\t" ;
    //            cout << endl;
    //            cout << " y =  " ;
    //            for (int i = 0; i < n; i++)
    //                cout << y[i] << "\t" ;
    //            cout << endl;
    //            cout << " just scale y, - - - - - - - - - - - - - current obj =  " << data.logDet_AplusX( y , data.OriginalL ) - data.logDet_AplusX( y , data.OriginalU )+ log(data.Sigma[n][n]);
    //            cout << endl;
    //            { int temp; cin >> temp; }
    //#endif
    //
    //
    //            // move y
    //            if (ite == ITSTEP)
    //                break;
    //            vector<double> grad = data.gradientLDYplusM(y, data.OriginalL );
    //            vector<double> gradU = data.gradientLDYplusM(y, data.OriginalU );
    //            vector<double> negativeShrinkedY(n);
    //            for (int i = 0; i < n; i++) {
    //                negativeShrinkedY[i] = -.95*y[i];
    //            }
    //#if test_adap_psd_cut_level2
    //            cout << " fawefrwagerhgerhrthutki" << endl;
    //#endif
    //            vector<double> gradCons = data.gradientLDYplusM(negativeShrinkedY, ConstraintLHS );
    //            double lambda = .01;
    //#if test_adap_psd_cut_level2
    //            cout << " fawefrwagerhgerhrthutki" << endl;
    //#endif
    //
    //
    //            for (int i = 0; i < n; i++)
    //                grad[i] = grad[i] - gradU[i] - lambda * gradCons[i] ;
    //#if test_adap_psd_cut_level2
    //            cout << " tag (rgegegeg4wtw4), just move y towards gradient direction " << endl;
    //            cout << " y =  " ;
    //            for (int i = 0; i < n; i++)
    //                cout << y[i] << "\t" ;
    //            cout << endl;
    //            cout << " grad =  " ;
    //            for (int i = 0; i < n; i++)
    //                cout << grad[i] << "\t" ;
    //            cout << endl;
    //            cout << " going to move y, - - -  ~  - - current obj =  " << data.logDet_AplusX( y , data.OriginalL ) - data.logDet_AplusX( y , data.OriginalU )+ log(data.Sigma[n][n]);
    //#endif
    //            for (int i = 0; i < n; i++)
    //                y[i] += STEPLENGTH * grad[i];
    //#if test_adap_psd_cut_level2
    //            cout << " new y =  " ;
    //            for (int i = 0; i < n; i++)
    //                cout << y[i] << "\t" ;
    //            cout << endl;
    //            cout << " just move y, - - -  ~  - - current obj =  " << data.logDet_AplusX( y , data.OriginalL ) - data.logDet_AplusX( y , data.OriginalU )+ log(data.Sigma[n][n]);
    //#endif
    //        }
    //
    //
    //
    //
    //#if test_adap_psd_cut_level2
    //        cout << " tag foaw9jf0aw34jf0aw39u0 ";
    //        cout << " y =  " ;
    //        for (int i = 0; i < n; i++)
    //            cout << y[i] << "\t" ;
    //        cout << endl;
    //        cout << " just move y, need to be scaled\n" ;
    //        { int temp; cin >> temp; }
    //#endif
    //
    //    }
    
#if test_adap_psd_cut_level3
    cout << "before scale up: " << endl;
    cout << " y = [ " ;
    for (int i = 0; i < n; i++)
        cout << y[i] << ", " ;
    cout << " ] " << endl;
#endif
    
    for (int i = 0; i < n; i++)
        sqrty[i] = sqrt(y[i]);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            LHS[i][j] = C[i][j] / sqrty[i] / sqrty[j];
    vector<double> e2 = data.eig(LHS) ;
    for (int i = 0; i < n; i++)
        y[i] *= e2[0];
    
#if test_adap_psd_cut_level3
    cout << "after scale up: " << endl;
    cout << " y = [ " ;
    for (int i = 0; i < n; i++)
        cout << y[i] << ", " ;
    cout << " ] " << endl;
    cout << " e2 = [ " ;
    for (int i = 0; i < n; i++)
        cout << e2[i] << ", " ;
    cout << " ] " << endl;
#endif
    
    
    
    
    vector<double> scale(n);
    for (int i = 0; i < n; i++) {
        scale[i] = sqrt( oneMinusXhatOverXhat[i] / y[i] );
    }
    
    
#if test_adap_psd_cut_level2
    
    cout << " tag (faw4q3b5sen) ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !  going to settle down the convex relaxation ! " << endl;
    
    
    cout << " faejfaewjfaew " << endl;
    cout << " data.OriginalInvOfInvSum = " << endl;
    for (int i = 0; i < n; i++){
        for (int j = 0; j< n; j++)
            cout << data.OriginalInvOfInvSum[i][j] << "\t" ;
        cout << endl;
    }
    
    cout << " y =  " ;
    for (int i = 0; i < n; i++)
        cout << y[i] << "\t" ;
    cout << endl;
    cout << " J =  " ;
    for (int i = 0; i < n; i++)
        cout << scale[i] << "\t" ;
    cout << endl;
#endif
    
    
    /*
     start to add cut !!!
     */
    IloRange cut;
    IloExpr rhs(env);
    try {
        vector < vector < double > > scaledL = data.scaleBothSide(data.OriginalL, scale);;
        vector < vector < double > > scaledU = data.scaleBothSide(data.OriginalU, scale);;
        
        
        
#if test_adap_psd_cut_level4
        cout << " \n\n\n\n tag (faw4zfg567tn65he4s) " << endl ;
        cout << " data.OriginalL = " << endl ;
        for (int i = 0; i < n; i++){
            for (int j = 0; j < n; j++)
                cout << data.OriginalL[i][j] << "\t" ;
            cout << endl;
        }
        cout << " data.OriginalU = " << endl ;
        for (int i = 0; i < n; i++){
            for (int j = 0; j < n; j++)
                cout << data.OriginalU[i][j] << "\t" ;
            cout << endl;
        }
        cout << " scaledL = " << endl ;
        for (int i = 0; i < n; i++){
            for (int j = 0; j < n; j++)
                cout << scaledL[i][j] << "\t" ;
            cout << endl;
        }
        cout << " scaledU = " << endl ;
        for (int i = 0; i < n; i++){
            for (int j = 0; j < n; j++)
                cout << scaledU[i][j] << "\t" ;
            cout << endl;
        }
        cout << " valz =  " ;
        for (int i = 0; i < n; i++)
            cout << valz[i] << "\t" ;
        cout << endl;
        cout << " scale =  " ;
        for (int i = 0; i < n; i++)
            cout << scale[i] << "\t" ;
        cout << endl;
        
        cout << " ready to add cut !! !! ! ! ! ! ! !  ! !! ! ! ! ~~ ~~~ ~ ~ ~ ~ ~ ~  ! ! ! !  ! !! ! ! ! ~~ ~~~ ~ ~ ~ ~ ~ ~  ! ! ! !  ! !! ! ! ! ~~ ~~~ ~ ~ ~ ~ ~ ~  ! ! ! !  ! !! ! ! ! ~~ ~~~ ~ ~ ~ ~ ~ ~  ! ! ! !  ! !! ! ! ! ~~ ~~~ ~ ~ ~ ~ ~ ~  ! ! ! !  ! !! ! ! ! ~~ ~~~ ~ ~ ~ ~ ~ ~  ! ! ! !  ! !! ! ! ! ~~ ~~~ ~ ~ ~ ~ ~ ~  ! ! ! !  ! !! ! ! ! ~~ ~~~ ~ ~ ~ ~ ~ ~  ! ! ! !  ! !! ! ! ! ~~ ~~~ ~ ~ ~ ~ ~ ~  ! ! ! !  ! !! ! ! ! ~~ ~~~ ~ ~ ~ ~ ~ ~  ! ! ! !  ! !! ! ! ! ~~ ~~~ ~ ~ ~ ~ ~ ~  ! ! ! !  ! !! ! ! ! ~~ ~~~ ~ ~ ~ ~ ~ ~  ! ! ! !  ! !! ! ! ! ~~ ~~~ ~ ~ ~ ~ ~ ~  ! ! ! !  ! !! ! ! ! ~~ ~~~ ~ ~ ~ ~ ~ ~  ! ! ! !  ! !! ! ! ! ~~ ~~~ ~ ~ ~ ~ ~ ~  ! ! ! !  ! !! ! ! ! ~~ ~~~ ~ ~ ~ ~ ~ ~  ! ! ! !  ! !! ! ! ! ~~ ~~~ ~ ~ ~ ~ ~ ~  =  " ;
        
        {
            int temp;
            cin >> temp;
            cout << " ! ! ! !  " ;
            cin >> temp;
        }
#endif
        
        
        
        double lpart = data.logDet_IplusDiagzW( valz, scaledL )  + log(data.Sigma[n][n]) ;
        double upart = data.logDet_IplusDiagzW( valz, scaledU ) ;
        double c = lpart - upart;
        
        /* test if current solution violates the cutting planes. if yes, add the cut */
        double criterior2 = (vall - valu - c);
#if test_adap_psd_cut
        cout << "   * adaptive psd cut * * adaptive psd cut * * adaptive psd cut * * adaptive psd cut * * adaptive psd cut * * adaptive psd cut * * adaptive psd cut * * adaptive psd cut * * adaptive psd cut *  " << endl;
        cout << "c = " << c << endl;
        cout << "lpart = " << lpart << ", upart = " << upart << ", lpart-upart = " << lpart-upart << endl;
        cout << "vall = " << vall << ", valu = " << valu << ", vall-valu = " << vall-valu << endl;
        cout << " ** criterior2 = (vall - valu - c) = " << criterior2 << endl;
        cout << " valz = "  << valz << endl;
#endif
        
        // rhs should be 0 now;
        if ( ( (nnode < 1000) && ( criterior2 < -1e-2 ) )  || ( (nnode > 1000) && ( criterior2 < -1e-1 ) ) ) {
            rhs += lpart - upart;
            vector<double> grad( data.uMinusLGradient( valz, scaledL, scaledU ) );
            for (int i = 0; i < n; i++){
                rhs += grad[i] * ( z[i] - valz[i] );
            }
            
            cut = ( l - u - rhs >= 0 );
#if test_adap_psd_cut
            cout << " cut = "  << cut << endl;
#endif
            
            cuts2BeAdded.add(cut);
            tune.num_gradient_cut_on_LMinusU_adaptive_psd++;
        }
        
        
        rhs.end();
        
    } catch (...) {
#if test_adap_psd_cut_level4
        cout << " \n\n\n\n tag (fsefwaegaer) " << endl ;
        cout << " ready to add cut !! !! ! ! ! ! ! !  ! !! ! ! ! ~~ ~~~ ~ ~ ~ ~ ~ ~  ! ! ! !  ! !! ! ! ! ~~ ~~~ ~ ~ ~ ~ ~ ~  ! ! ! !  ! !! ! ! ! ~~ ~~~ ~ ~ ~ ~ ~ ~  ! ! ! !  ! !! ! ! ! ~~ ~~~ ~ ~ ~ ~ ~ ~  ! ! ! !  ! !! ! ! ! ~~ ~~~ ~ ~ ~ ~ ~ ~  ! ! ! !  ! !! ! ! ! ~~ ~~~ ~ ~ ~ ~ ~ ~  ! ! ! !  ! !! ! ! ! ~~ ~~~ ~ ~ ~ ~ ~ ~  ! ! ! !  ! !! ! ! ! ~~ ~~~ ~ ~ ~ ~ ~ ~  ! ! ! !  ! !! ! ! ! ~~ ~~~ ~ ~ ~ ~ ~ ~  ! ! ! !  ! !! ! ! ! ~~ ~~~ ~ ~ ~ ~ ~ ~  ! ! ! !  ! !! ! ! ! ~~ ~~~ ~ ~ ~ ~ ~ ~  ! ! ! !  ! !! ! ! ! ~~ ~~~ ~ ~ ~ ~ ! ! ! ! ! ! ! ~~ ~~~ ~ ~ ~ ~ ~ ~  ! ! ! !  ! !! ! ! ! ~~ ~~~ ~ ~ ~ ~ ~ ~  ! ! ! !  ! !! ! ! ! ~~ ~~~ ~ ~ ~ ~ ~ ~  ! ! ! !  ! !! ! ! ! ~~ ~~~ ~ ~ ~ ~ ~ ~  ! ! ! !  ! !! ! ! ! ~~ ~~~ ~ ~ ~ ~ ~ ~  =  " ;
        
        {
            int temp;
            cout << " ! ! ! !  " ;
            cin >> temp;
        }
#endif
        rhs.end();
        throw;
    };
    
};























//            ########     ###    ########  ######## ####    ###    ##          ########      ######   ########     ###    ########
//            ##     ##   ## ##   ##     ##    ##     ##    ## ##   ##          ##     ##    ##    ##  ##     ##   ## ##   ##     ##
//            ##     ##  ##   ##  ##     ##    ##     ##   ##   ##  ##          ##     ##    ##        ##     ##  ##   ##  ##     ##
//            ########  ##     ## ########     ##     ##  ##     ## ##          ########     ##   #### ########  ##     ## ##     ##
//            ##        ######### ##   ##      ##     ##  ######### ##          ##   ##      ##    ##  ##   ##   ######### ##     ##
//            ##        ##     ## ##    ##     ##     ##  ##     ## ##          ##    ##     ##    ##  ##    ##  ##     ## ##     ##
//            ##        ##     ## ##     ##    ##    #### ##     ## ########    ##     ##     ######   ##     ## ##     ## ########
void SparseRegressionProblem :: generateCuts_partial_r_gradient( IloEnv &env, IloConstraintArray &cuts2BeAdded, DataModel &data, TuneParameters &tune, IloBoolVarArray &z, IloNumVar &u, IloNumVar &l, IloNumArray valz, IloNum valu, IloNum vall, long nnode, vector<double>  &lastValz)
{
    
    int n = data.getCols();
    double r = data.rStar;
    //    if ( (tune.num_gradient_upper_cut_on_LMinusU > n * 3 ) && (nnode == 0) )
    //        return;
    //    if ( (nnode > 1000) && (rand()%10 > 0 ) )
    //        return;
    
    double fractional_part = 0;
    for (int i = 0; i < n; i++){
        fractional_part += min( abs(valz[i]), abs(valz[i]-1) );
    }
    if (fractional_part < 1e-4 * n)
        return;
    
    
    
    double distFromLastValz = 0;
    for (int i = 0; i < n ; i++){
        ;
        distFromLastValz += pow(lastValz[i] - valz[i], 2);
        
    }
    cout << " distFromLastValz = " << distFromLastValz << endl;
    if (distFromLastValz < .01*n)
        return;
    
    double epsilon = .0;
    
    IloNumArray realValz(env,n);
    for (int i = 0; i < n; i++){
        realValz[i] = valz[i];
        
        if (valz[i] < epsilon )
            valz[i] = epsilon;
        else if (valz[i] > 1.0 - epsilon)
            valz[i] = 1.0 - epsilon;
    }
    
    IloRange cut;
    IloExpr rhs(env);
    // u <= logdet(I + W * diag(z)) , W = Sigma - I
    // u <= logdet(I + W * diag(zhat))  + gradient' * (z - zhat)
    try {
        double lpart = r*r *  ( data.logDet_IplusDiagzW( valz, data.L )  + log(data.Sigma[n][n]) );
        double upart = data.logDet_IplusDiagzW( valz, data.U ) ;
        double c = lpart - upart;
        //        vector<double> grad( data.uMinusLGradient( valz ) );
        vector<double> grad( data.R2uMinusLGradient( valz, r, data.LMaxEigen1, data.UMaxEigen1 ) );
        for (int i = 0; i < n; i++){
            c += grad[i] * (realValz[i] - valz[i]);
        }
        
        /* test if current solution violates the cutting planes. if yes, add the cut */
        double criterior = ((r*r * vall - valu - c)/max(vall, valu));
        double criterior2 = (r*r * vall - valu - c);
        
        //        if ( criterior2 < - .01 ) {
        if ( criterior2 < -1e-2 ) {
            rhs += r*r*lpart - upart;
            for (int i = 0; i < n; i++){
                rhs += grad[i] * ( z[i] - valz[i] );
            }
            cut = ( r*r*l - u - rhs >= 0 );
            //        if ( criterior < -.0001 ) {
            cuts2BeAdded.add(cut);
            tune.num_gradient_upper_cut_on_LMinusU_partial++;
        }
        
        
#if pauseUpperGrad
        cout << " * const = " << c << endl;
        cout << " * grad = [ ";
        for (int i = 0; i < n; i++)
            cout << grad[i] <<", " ;
        cout << " ]\n";
        
        cout << " * valz = ";
        cout << valz ;
        cout << "\n";
        
        cout << " * cut = { "<< cut << "} "  << endl;
        cout <<" - - (upper gradient) {tag f34fq34f} - - -  "<< endl;
        if ( 1 ){
            int temp;
            cin >> temp;
        }
#endif
        rhs.end();
    } catch (...) {
        rhs.end();
        throw;
    }
};








