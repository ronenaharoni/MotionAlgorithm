#include <math.h>
#include <stdlib.h>
#include "MotionHeader.h"
#include "polyfit.h"

int MotionTracking(float* Mscan[],FgParams_Struct* FgParams,BgRadarParams *bgParams){
	int i,k,t,m;
	int DistThr= 5;
	float Centoird1,Centoird2, Noise_dB;
	float SumOfSecondDiff,MeanOfDiff,FirstDiff[FgParams->Motion.Nscans-1],SecondDiff[FgParams->Motion.Nscans-2];
	float SumForVar,NoiseEstimation,PowerOfMscan,maxEnergyPerTime;
	//	float energyPerTime_dB[FgParams->Motion.Nbins][FgParams->Motion.SpectrogramTimeBinsSingleMotion];
	double 	ValuesForInterp[FgParams->Motion.SpectrogramTimeBinsSingleMotion];
	double PointsForInterp[FgParams->Motion.SpectrogramTimeBinsSingleMotion];
	int NumPointsForInterp=0;
	double p[2];//1st order polynom coefficients
	float   NextRange,RstartWanted;
	float* energyPerTime_dB[FgParams->Motion.Nbins];//,*Mscan_old[FgParams.Nscans];
	for (i=0; i<FgParams->Motion.Nbins; i++){
		energyPerTime_dB[i] = (float *)malloc(FgParams->Motion.SpectrogramTimeBinsSingleMotion * sizeof(float));
	}

	//	int FirstBinAboveThresholdPerTime[FgParams->Motion.SpectrogramTimeBinsSingleMotion];
	//	memset(FirstBinAboveThresholdPerTime,0,(FgParams->Motion.SpectrogramTimeBinsSingleMotion)*sizeof(int));

	for(i=0;i<FgParams->Motion.Nbins;i++){//calculate NoiseEstimation=var(diff(Signal_k,2)) for each bin
		SumOfSecondDiff=0;
		for(k=0;k<FgParams->Motion.Nscans-1;k++){//1st Diff
			FirstDiff[k]=(Mscan[k+1][i]-Mscan[k][i]);
		}
		for(k=0;k<FgParams->Motion.Nscans-2;k++){//2nd Diff
			SecondDiff[k]=(FirstDiff[k+1]-FirstDiff[k]);
			SumOfSecondDiff+=SecondDiff[k];
			//			printf("%d %lf\n",k,SecondDiff[k]);
		}

		MeanOfDiff=SumOfSecondDiff/(FgParams->Motion.Nscans-2);//Calculate mean
		SumForVar=0;
		for(k=0;k<FgParams->Motion.Nscans-2;k++){//VAR=SUM((x-mu)^2)/N-1,  N-1 because it's unbiased estimator
			SumForVar+=pow((SecondDiff[k]-MeanOfDiff),2);
		}
		NoiseEstimation=SumForVar/((FgParams->Motion.Nscans-2)-1);//it's the Variance

		for(t=0;t<FgParams->Motion.SpectrogramTimeBinsSingleMotion;t++){//take 30 points of slow and calculate the max energy of them
			maxEnergyPerTime=0;
			for(m=0;m<FgParams->Motion.SpectrogramWinLength;m++){
				PowerOfMscan=pow(Mscan[m+t*(FgParams->Motion.SpectrogramTimeShift)][i],2);
				//				printf("%d %lf\n",m,Mscan[m+t*(FgParams->Motion.TimeShift)][i]);

				if(PowerOfMscan > maxEnergyPerTime){
					maxEnergyPerTime=PowerOfMscan;
				}
			}

			energyPerTime_dB[i][t]=10*log10(maxEnergyPerTime/NoiseEstimation);//maybe do minus instead/
			//			printf("%d %d    %lf\n",i,t,energyPerTime_dB[i][t]);

		}

	}

	K_means(energyPerTime_dB,FgParams,&Centoird1,&Centoird2);//Calulcate the centroids of the samples

	if(fabs(Centoird1-Centoird2)>DistThr){
		if(Centoird1>Centoird2)
			Noise_dB =Centoird1;
		else
			Noise_dB =Centoird2;
	}
	else
		Noise_dB=INFINITY;//why??



	//1st order interpolation for the curve of the walking

	for(t=0;t<FgParams->Motion.SpectrogramTimeBinsSingleMotion;t++){
		for(i=0;i<FgParams->Motion.Nbins;i++){
			if(energyPerTime_dB[i][t]>Noise_dB){//find the first Bin which crosses the energy
				ValuesForInterp[NumPointsForInterp]=FgParams->Motion.Rbin_m[i];//this is the first rbin_m who cross the threshold put in rbin_m
				PointsForInterp[NumPointsForInterp]=t+1;//save the valid time point for interpolation, +1 to be like matlab
				NumPointsForInterp+=1;
				break;
			}
		}

	}

	//Add this
	if((FgParams->Motion.SpectrogramTimeBinsSingleMotion - NumPointsForInterp) > 0.87*FgParams->Motion.SpectrogramTimeBinsSingleMotion){//if we don't have enough points -> don't change
//		FgParams->Rstart_corrected=FgParams->Rbin_m[0];
//		FgParams->Rstop_corrected=FgParams->Rbin_m[FgParams->Motion.SpectrogramTimeBinsSingleMotion-1]-0.1;
		bgParams->Rstart_m=FgParams->Motion.Rbin_m[0];
		bgParams->Rstop_m=FgParams->Motion.Rbin_m[FgParams->Motion.SpectrogramTimeBinsSingleMotion-1]-0.1;


			for (i = 0; i < FgParams->Motion.Nbins; i++) {
				free(energyPerTime_dB[i]);
			}

		return 0;
	}

	polyfit(PointsForInterp,ValuesForInterp,NumPointsForInterp,1,p);//1 si the order MUST TO BE float

	NextRange=p[1]*FgParams->Motion.SpectrogramTimeBinsSingleMotion+p[0];//insert the last point in order to detemine the center of the next ring

	RstartWanted=bgParams->Rmin;
	if ((NextRange-0.5*2.5) > RstartWanted)//max(nextRange - 1.25, Rmin)
		RstartWanted=NextRange-0.5*2.5;

	if((RstartWanted + 2.5) > bgParams->Rmax){//if the end of the ring is exit from the Rmax, cut it by the Rmax
		bgParams->Rstart_m = bgParams->Rmax - 2.5;
		bgParams->Rstop_m=bgParams->Rmax;
	}
	else{
		bgParams->Rstart_m = RstartWanted;
		bgParams->Rstop_m = RstartWanted + 2.5;
	}

	for (i = 0; i < FgParams->Motion.Nbins; i++) {
				free(energyPerTime_dB[i]);
			}

	return 0;
}


int K_means(float* energyPerTime_dB[],FgParams_Struct* FgParams,float* ChosenCentoird1,float* ChosenCentoird2){
	float SumForCenteroid1, SumForCenteroid2,Centeroid1,Centeroid2,PrevCenteroid1,PrevCenteroid2;
	 float epsilon=1e-5;
	int i,t,AlgorithmDone=0;
	int NumOfSamplesCentroid1,NumOfSamplesCentroid2,NumOfIter=0;
	SumForCenteroid1=0;
	SumForCenteroid2=0;
	int Half_K_mean_TotalSamples=FgParams->Motion.SpectrogramTimeBinsSingleMotion * FgParams->Motion.Nbins/2;
	//find initial centroids by deviding to two groups
	for(i=0;i<FgParams->Motion.Nbins;i++){
		for(t=0;t<FgParams->Motion.SpectrogramTimeBinsSingleMotion;t++){
			if(t<FgParams->Motion.SpectrogramTimeBinsSingleMotion/2)//first half
				SumForCenteroid1+=energyPerTime_dB[i][t];
			else//second half
				SumForCenteroid2+=energyPerTime_dB[i][t];
		}
	}
	Centeroid1=SumForCenteroid1/(Half_K_mean_TotalSamples);
	Centeroid2=SumForCenteroid2/(Half_K_mean_TotalSamples);


	while(AlgorithmDone==0){
		SumForCenteroid1=0;
		SumForCenteroid2=0;
		NumOfSamplesCentroid1=0;
		NumOfSamplesCentroid2=0;

		for(i=0;i<FgParams->Motion.Nbins;i++){
			for(t=0;t<FgParams->Motion.SpectrogramTimeBinsSingleMotion;t++){

//				printf("%d %d    %lf\n",i,t,energyPerTime_dB[i][t]);

				if(fabs(energyPerTime_dB[i][t] - Centeroid1) < fabs(energyPerTime_dB[i][t] - Centeroid2)){//if the sample is closer to centroid 1
					SumForCenteroid1+=energyPerTime_dB[i][t];
					NumOfSamplesCentroid1+=1;
				}
				else{//the sample is closer to centroid2
					SumForCenteroid2+=energyPerTime_dB[i][t];
					NumOfSamplesCentroid2+=1;
				}
			}
		}
			//new centroids:
			PrevCenteroid1=Centeroid1;
			PrevCenteroid2=Centeroid2;

			Centeroid1=SumForCenteroid1/NumOfSamplesCentroid1;
			Centeroid2=SumForCenteroid2/NumOfSamplesCentroid2;


		//Stop condition
		if ((fabs(Centeroid1-PrevCenteroid1) <= epsilon) && (fabs(Centeroid2-PrevCenteroid2) <= epsilon)){
			*ChosenCentoird1=Centeroid1;
			*ChosenCentoird2=Centeroid2;
			AlgorithmDone=1;
		}
		NumOfIter++;
	}
	return 0;

}
