#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <ctime>
#include <math.h>
#include <cstring>
#include <random>
#include <chrono>
#include <functional>
#include <map>
#include <vector>
#include <array>
#include <algorithm>
#include <functional>

using namespace std;

//TUNABLE PARAMETERS
//Base
const string version {"Version4.x Generic"};
const int SPACESIZE {60'000'000};
const int MAXDAYS  {180};
const int START {10'000};
const int OUTPUT_TYPE{0};
const double SAMPLE_FRACTION {0.01};
const long RSEED {12};
//BrainGeometry
const int SEEDING_SIZE {157'000};
const int BASE_VOXEL {50'000};
const double SCALE_VOX {6};
int brainVoxel {BASE_VOXEL*SCALE_VOX};
const int QUBI {200};
//Cell Cycles parameters are given as hours*10^(-2)
const int LOW_CYCLE_LIMIT {850};
const int BASE_CYCLE {-2200};
const double SD_BASE_CYCLE {-80};
const int SD_DELTA_CYCLE {250};
//skewnorm parameters from timelaspe data
const double SEEDCYCLE_MODE {13.5};
const double SEEDCYCLE_SCALE{9};
const double SEEDCYCLE_SHAPE{11};
//Migration
const int MIGRATION_UNIT {7500};
//CellDeath
const double P_DEATH_SCALE {0.2};
const double EC50 {60};
const double P_DEATH_T_ZERO {0.0025};
const double ALIEN_WEIGHT {1};
const double LINEAGE_WEIGHT {1};
//Quiescence and Stemness
const double G_ZERO_P {0.45};
const double STEM_FREQ_T0 {1};
const int N_DIFF_CYCLE {12};
const int SENCYCL{2};
const double P_ASYM {0}; //probability of asymmetric cell division
const double P_T_GEN {0};//probability of no two stem cells
//END TUNABLE PARAMETERS



//NON-TUNABLE GLOBALS
uint_fast64_t mseed{0};
//END OF NON-TUNABLE GLOBALS




//SAVE PARAMETER SUMMARY IN A FILE "XXX_1-PARAMETERS.TXT"
int printParameter(std::string  ParamFileName) {

    std::ofstream fparam;
    fparam.open (ParamFileName);

    if(!fparam.is_open()) return -1;

    fparam<<"FILE NAME: "<<ParamFileName<<endl;
    fparam<<"Version: "<<version<<endl;
    fparam<<"***General*** "<<endl;
    fparam<<"SPACESIZE "<<SPACESIZE<<" cells"<<endl;
    fparam<<"MAXDAYS "<<MAXDAYS<<" days"<<endl;
    fparam<<"START "<<START<<" cells"<<endl;
    fparam<<"OUTPUT_TYPE "<<OUTPUT_TYPE<<" cells"<<endl;
    fparam<<"SAMPLE_FRACTION "<<SAMPLE_FRACTION<<" cells"<<endl;
    fparam<<"***BrainGeometry***"<<endl;
    fparam<<"SEEDING_SIZE "<<SEEDING_SIZE<<" nm"<<endl;
    fparam<<"brainVoxel "<<brainVoxel<<" nm"<<endl;
    fparam<<"SCALE_VOX "<<SCALE_VOX<<" (RATIO)"<<endl;
    fparam<<"QUBI "<<QUBI<<endl;
    fparam<<"***CellCycle***"<<endl;
    fparam<<"LOW_CYCLE_LIMIT "<<LOW_CYCLE_LIMIT<<" centiHrs"<<endl;
    fparam<<"BASE_CYCLE "<<BASE_CYCLE<<" centiHrs"<<endl;
    fparam<<"SD_BASE_CYCLE "<<SD_BASE_CYCLE<<endl;
    fparam<<"SD_DELTA_CYCLE "<<SD_DELTA_CYCLE<<endl;
    fparam<<"SEEDCYCLE_MODE "<<SEEDCYCLE_MODE<<endl;
    fparam<<"SEEDCYCLE_SCALE "<<SEEDCYCLE_SCALE<<endl;
    fparam<<"SEEDCYCLE_SHAPE "<<SEEDCYCLE_SHAPE<<endl;
    fparam<<"***Migration***"<<endl;
    fparam<<"MIGRATION_UNIT "<<MIGRATION_UNIT<<" nm"<<endl;
    fparam<<"***CellDeath***"<<endl;
    fparam<<"P_DEATH_SCALE "<<P_DEATH_SCALE<<endl;
    fparam<<"EC50 "<<EC50<<endl;
    fparam<<"P_DEATH_T_ZERO "<<P_DEATH_T_ZERO<<endl;
    fparam<<"ALIEN_WEIGHT "<< ALIEN_WEIGHT<<endl;
    fparam<<"LINEAGE_WEIGHT "<<LINEAGE_WEIGHT<<endl;
    fparam<<"***Quiescence and stemness***"<<endl;
    fparam<<"G_ZERO_P "<<G_ZERO_P<<endl;
    fparam<<"P_ASYM  "<<P_ASYM <<endl;
    fparam<<"P_T_GEN "<<P_T_GEN<<endl;
    fparam<<"SENCYCL "<<SENCYCL<<endl;
    fparam<<"RSEED "<<RSEED<<endl;
    fparam<<"mSEED "<<mseed<<endl;
    fparam.close();
    return 1;
}
//END SAVE PARAMETER SUMMARY

//PSEUDORANDOM GENERATOR FOR SKEWNORM DISTRIBUTION
auto cfnorm= std::bind(std::normal_distribution<double>(0,1),mt19937_64 (mseed));

double rskn(double loc, double scale, double shape) {

    double rho = shape/sqrt(1 + shape*shape);
    double u0 =cfnorm();
    double v = cfnorm();
    double u1=rho * u0 + sqrt(1 - rho*rho) * v;
    return(10*(loc + (scale * u1)*(u0/abs(u0))));
}
//END PSEUDORANDOM GENERATOR FOR SKEWNORM DISTRIBUTION


int main(int argc, char* argv[]) {

    std::hash<long> hashlong;
    mseed= hashlong(-RSEED);

    auto cunif = std::bind(std::uniform_real_distribution<double>(0.0,1.0),mt19937_64 (mseed));
    auto cnorm= std::bind(std::normal_distribution<double>(0,1),mt19937_64 (mseed));

    std::string TrackByHrsFileName;
    std::string FinalPopFileName;
    std::string CloneTableFileName;



    auto t = std::time(nullptr);
    auto tm = *std::localtime(&t);
    std::ostringstream oss;
    oss << std::put_time(&tm, "%d_%m_h%H_%M_");


    printParameter(oss.str()+"_1-Parameters.txt");


    TrackByHrsFileName=oss.str()+"_2-TrackByHours.txt";
    std::ofstream ftrackH;
    ftrackH.open (TrackByHrsFileName);
    ftrackH<<"days"<<"\t";
    ftrackH<<"popsize"<<"\t";
    ftrackH<<"nclones"<<"\t";
    ftrackH<<"shannon"<<"\t";
    ftrackH<<"CyleMedian"<<"\t";
    ftrackH<<"restingCells"<<"\t";
    ftrackH<<"numberOfDeaths"<<"\t";
    ftrackH<<"newCells"<<"\t";
    ftrackH<<"lenProliferating"<<"\t";
    ftrackH<<"stemNumber"<<"\t";
    ftrackH<<"time"<<endl;

    if(!ftrackH.is_open()) return -1;

    FinalPopFileName=oss.str()+"_4-FinalPopulation.txt";

    CloneTableFileName=oss.str()+"_3-CloneTableBy10dd.txt";
    std::ofstream fclonetable;
    fclonetable.open (CloneTableFileName);


    int debugg=0;

    std::ofstream debugtable;

	//DEBUGTABLE FILE CREATION
    if(debugg) {
        std::string debugTableFileName;
        debugTableFileName=oss.str()+"_5-debugTableDensity.txt";
        std::ofstream debugtable;
        debugtable.open (debugTableFileName);
    }
	//END DEBUGTABLE FILE CREATION



	//GENERATE START TIME STRING
    t = std::time(nullptr);
    tm = *std::localtime(&t);
    oss.str("");
    oss << std::put_time(&tm, "@ %H:%M:%S");
	//END GENERATE START TIME STRING



    const int quby{QUBI};
    const int qubz{quby*quby};
    const int sizeDensM{quby*quby*quby};
    std::vector <map <int, double>> densM(sizeDensM);
    cout<<"size densM: "<<densM.size()<<endl;

    vector<array <int,7>> cellsp(SPACESIZE, {0,0,0,0,0,0,0});
    vector<array <int,7>> cellsd(SPACESIZE, {0,0,0,0,0,0,0});

    int nclones=0;
    int popsize=0;
    int numberOfDeaths=0;
    int newCells=0;
    double CyleMedian=0;
    double shannon=0;
    int restingCells=0;
    map <int,int> cloneMap;
    vector <double> cycles;
    int stemNumber{0};



	//CELL SEEDING LOOP
    for(int i=0; i<START; i++) {
        cellsp[i][0]=i+1;
        cellsp[i][2]=10*rskn(SEEDCYCLE_MODE,SEEDCYCLE_SCALE,SEEDCYCLE_SHAPE);
        cellsp[i][1]=cunif()*cellsp[i][2];
        cellsp[i][3]=(cunif()* 2*SEEDING_SIZE)-SEEDING_SIZE;//x
        cellsp[i][4]=(cunif()* 2*SEEDING_SIZE)-SEEDING_SIZE;//y
        cellsp[i][5]=(cunif()* 4*SEEDING_SIZE)-2*SEEDING_SIZE;//z
        cellsp[i][6]= (cunif()<STEM_FREQ_T0)? (1):(-SENCYCL);
    }
	//END CELL SEEDING LOOP



	//START SIMULATION LOOP

    popsize=START;
    int xi,yi,zi;

        for(int simulatedDays=0; simulatedDays<MAXDAYS; simulatedDays++) {
            for(int simulatedHours=0; simulatedHours<24; simulatedHours++) {


				//CALCULATING CELL DENSITY MAP
                for(int i=0; i<popsize; i++) {
                    xi=(cellsp[i][3]+brainVoxel*(QUBI/2))/brainVoxel;
                    yi=(cellsp[i][4]+brainVoxel*(QUBI/2))/brainVoxel;
                    zi=(cellsp[i][5]+brainVoxel*(QUBI/2))/brainVoxel;
                    xi=(xi>(QUBI-1))? (QUBI-1):xi;
                    xi=(xi<0)? 0:xi;
                    yi=(yi>(QUBI-1))? (QUBI-1):yi;
                    yi=(yi<0)? 0:yi;
                    zi=(zi>(QUBI-1))? (QUBI-1):zi;
                    zi=(zi<0)? 0:zi;
                    int densMcoord=xi+quby*yi+qubz*zi;
                    densM.at(densMcoord)[cellsp[i][0]]++;
                    densM.at(densMcoord)[0]++;
                }
				//END CALCULATING CELL DENSITY MAP


                int Ndcell=0;

				//SURVIVORS IDENTIFICATION AND COPY
                for(int pcell=0; pcell<popsize; pcell++) {

                    xi=(cellsp[pcell][3]+brainVoxel*(QUBI/2))/brainVoxel;
                    yi=(cellsp[pcell][4]+brainVoxel*(QUBI/2))/brainVoxel;
                    zi=(cellsp[pcell][5]+brainVoxel*(QUBI/2))/brainVoxel;
                    xi=(xi>(QUBI-1))? (QUBI-1):xi;
                    xi=(xi<0)? 0:xi;
                    yi=(yi>(QUBI-1))? (QUBI-1):yi;
                    yi=(yi<0)? 0:yi;
                    zi=(zi>(QUBI-1))? (QUBI-1):zi;
                    zi=(zi<0)? 0:zi;
                    int densMcoord=xi+quby*yi+qubz*zi;
                    int alienDensity=(densM.at(densMcoord)[0])-(densM.at(densMcoord)[cellsp[pcell][0]]);
                    int autoDensity=(densM.at(densMcoord)[cellsp[pcell][0]]);

                    double pDeath=P_DEATH_T_ZERO+(1/(1+exp(-P_DEATH_SCALE
                                                           *((1/(SCALE_VOX*SCALE_VOX*SCALE_VOX))*(alienDensity*ALIEN_WEIGHT+autoDensity*LINEAGE_WEIGHT)-EC50)))
                                                  *(1-P_DEATH_T_ZERO));


                    if(cunif()>pDeath) {
                        cellsd[Ndcell][0]=cellsp[pcell][0];
                        cellsd[Ndcell][1]=cellsp[pcell][1]-100;
                        cellsd[Ndcell][2]=cellsp[pcell][2];
                        cellsd[Ndcell][3]=cellsp[pcell][3]+cnorm()*MIGRATION_UNIT;//HOURLY MOVEMENT - X
                        cellsd[Ndcell][4]=cellsp[pcell][4]+cnorm()*MIGRATION_UNIT;//HOURLY MOVEMENT - Y
                        cellsd[Ndcell][5]=cellsp[pcell][5]+cnorm()*MIGRATION_UNIT;//HOURLY MOVEMENT - Z
                        cellsd[Ndcell][6]=cellsp[pcell][6];
                        Ndcell++;
                    }
                } 
				// END SURVIVORS IDENTIFICATION AND COPY

                numberOfDeaths=popsize-Ndcell;
                newCells=0;
                restingCells=0;

				//REPLICATION LOOP

                for (int dcell=0; dcell<Ndcell; dcell++) {

                    int newcellIndex=Ndcell+newCells;

                    if(cellsd[dcell][1]<=0) {
                        if(newcellIndex>SPACESIZE) {
                            ftrackH.close();
                            fclonetable.close();
                            return -1;
                        }

                        cellsd[newcellIndex][0]=cellsd[dcell][0];
						
						//CELL CYCLE GENERATION FOR DAUGHTER CELL 2
                        double gauss=cnorm()*SD_DELTA_CYCLE;
                        int goG0=(cunif()<=G_ZERO_P)||(cellsd[dcell][6]==0);
                        double newCycleLen=cellsd[dcell][2]             
                                           +(gauss<0)*gauss*((cellsd[dcell][2]/LOW_CYCLE_LIMIT)-1)
                                           +(gauss>=0)*gauss+1e7*goG0;
						//END CELL CYCLE GENERATION FOR DAUGHTER CELL 2
						
                        cellsd[newcellIndex][1]= newCycleLen;
                        cellsd[newcellIndex][2]= newCycleLen;

                        cellsd[newcellIndex][3]=cellsd[dcell][3]+cnorm()
                                                *(MIGRATION_UNIT/5); //HOURLY MOVEMENT OF DAUGHTER CELL - X

                        cellsd[newcellIndex][4]=cellsd[dcell][4]+cnorm()
                                                *(MIGRATION_UNIT/5); //HOURLY MOVEMENT OF DAUGHTER CELL - Y

                        cellsd[newcellIndex][5]=cellsd[dcell][5]+cnorm()
                                                *(MIGRATION_UNIT/5); //HOURLY MOVEMENT OF DAUGHTER CELL - Z


                        cellsd[newcellIndex][6]=cellsd[dcell][6]+1;

						//CELL CYCLE GENERATION FOR DAUGHTER CELL 1
                        gauss=cnorm()*SD_DELTA_CYCLE;
                        goG0=(cunif()<=G_ZERO_P)||(cellsd[dcell][6]==0);
                        restingCells+=goG0;
                        newCycleLen=cellsd[dcell][2]
                                    +(gauss<0)*gauss*((cellsd[dcell][2]/LOW_CYCLE_LIMIT)-1)
                                    + (gauss>=0)*gauss+1e7*goG0;
						//END CELL CYCLE GENERATION FOR DAUGHTER CELL 1
						
                        cellsd[dcell][2]=cellsd[dcell][1]= newCycleLen;
                        cellsd[dcell][6]=goG0?(-100):(cellsd[dcell][6]+1);

                        if(cellsd[dcell][6]>0) {
                            double divisionType=cunif();

                            if (divisionType<=P_ASYM) {
                                cellsd[newcellIndex][6]=-N_DIFF_CYCLE;

                            } else if(divisionType<=P_T_GEN) {
                                cellsd[newcellIndex][6]=-N_DIFF_CYCLE;
                                cellsd[dcell][6]=-N_DIFF_CYCLE;
                            }
                        }

                        newCells++;
                    }
                }
				//END REPLICATION LOOP

                int maxsize=(((popsize)>(newCells+Ndcell))?(popsize):(newCells+Ndcell));
                popsize=newCells+Ndcell;

                //////////////////// COPY FROM CELLSD TO CELLSP
                for (int j=0; j<7; j++) {
                    for(int i=0; i<maxsize; i++) {
                        cellsp[i][j]=cellsd[i][j];
                    }
                }
                 ////////////////////END COPY FROM CELLSD TO CELLSP

                cloneMap.clear();
                cycles.clear();
                stemNumber=0;

                for(int i=0; i<popsize; i++) {
                    cloneMap[cellsp[i][0]]++;

                    if(cellsp[i][2]<1e5) cycles.push_back(cellsp[i][2]);

                    if(cellsp[i][6]>0) stemNumber++;
                }

                int lenProliferating=cycles.size();

                shannon=0;
                nclones=cloneMap.size();

                for(auto clone:cloneMap) {

                    shannon=shannon+((double)clone.second/(double)popsize)
                            *log((double)(clone.second)/(double)(popsize));
                }

                shannon=-shannon;

                int n;
                n=cycles.size()/2;
                nth_element(cycles.begin(),cycles.begin()+n,cycles.end());
                CyleMedian=cycles[n];

                densM.clear();
                densM.resize(sizeDensM);

                t = std::time(nullptr);
                tm = *std::localtime(&t);
                oss.str("");
                oss << std::put_time(&tm, "@ %H:%M:%S");

                //HOUR-BY-HOUR TRACKING
                ftrackH<<simulatedDays<<"\t";
                ftrackH<<popsize<<"\t";
                ftrackH<<nclones<<"\t";
                ftrackH<<shannon<<"\t";
                ftrackH<<CyleMedian/100<<"\t";
                ftrackH<<restingCells<<"\t";
                ftrackH<<numberOfDeaths<<"\t";
                ftrackH<<newCells<<"\t";
                ftrackH<<lenProliferating<<"\t";
                ftrackH<<stemNumber<<"\t";
                ftrackH<<oss.str()<<endl;
				//END HOUR-BY-HOUR TRACKING
				
            } //HOUR-LOOP END

            //DAY-BY-DAY TRACKING IN CONSOLE
            cout<<simulatedDays<<"\t";
            cout<<popsize<<"\t";
            cout<<nclones<<"\t";
            cout<<CyleMedian/100<<"\t";
            cout<<stemNumber<<"\t";
            cout<<oss.str()<<endl;
            //END DAY-BY-DAY TRACKING IN CONSOLE

			//CLONE TABLE RECORDING EACH 10 DAYS
            if(simulatedDays%10==0) {
                for(auto i :cloneMap)    {
                    fclonetable<<i.first<<"\t"<<i.second<<endl;
                }

                fclonetable<<endl;
                fclonetable.flush();

                std::ofstream fFinalPop;
                fFinalPop.open (FinalPopFileName);

                if (fFinalPop.is_open()) {

                    for(int i=0; i<popsize; i++) {
                        if(cunif()<SAMPLE_FRACTION||OUTPUT_TYPE) {
                            fFinalPop<<cellsp[i][0]<<"\t";
                            fFinalPop<<cellsp[i][1]<<"\t";
                            fFinalPop<<cellsp[i][2]<<"\t";
                            fFinalPop<<cellsp[i][3]<<"\t";
                            fFinalPop<<cellsp[i][4]<<"\t";
                            fFinalPop<<cellsp[i][5]<<"\t";
                            fFinalPop<<cellsp[i][6]<<"\n";
                        }
                    }

                    fFinalPop.close();

                } else {
                    printf("Failed to create the files.\n");
                    return -1;
                }


            }
			//END CLONE TABLE RECORDING EACH 10 DAYS
			
			
            ftrackH.flush();
        }
		
		/////END SIMULATION LOOP//////////////



    ftrackH.close();
    fclonetable.close();



    return 0;
}
