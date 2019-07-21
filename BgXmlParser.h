/*
 * BgXmlParser.h
 *
 *  Created on: Mar 7, 2019
 *      Author: echocare
 */

#ifndef BGXMLPARSER_H_
#define BGXMLPARSER_H_
//#include "hostInterfaceCommon.h"

#ifdef __cplusplus
extern "C" {
#endif






typedef struct {
	int Nscans;
	int Nbins;
	int Fs;
	int RF_NumOfClasses;
	int SpectrogramRangeWidth;
	int SpectrogramWinLength;
	int SpectrogramN_Overlap;
	int SpectrogramNoiseFreqBins;
	int FminBin;
	int Fbins;
	int NoiseThresh;
	int TopMaxFreq;
	int MinFreqPlusMinus;
	float MinFreqReal;
	int MinEventDuration;
	int GapLength;
	int NumSamplesForDerivativeEstimation;
	int MedianValue;
	int AvgValue;
	int TruncateReal;
	int TruncateHilbert;

} MotionXML_Params;




typedef enum
{
	Nscans                  ,
	Nbins                   ,
	Fs                      ,
	RF_NumOfClasses         ,
	SpectrogramRangeWidth   ,
	SpectrogramWinLength    ,
	SpectrogramN_Overlap    ,
	SpectrogramNoiseFreqBins,
	FminBin                 ,
	Fbins                   ,
	NoiseThresh             ,
	TopMaxFreq              ,
	MinFreqPlusMinus        ,
	MinFreqReal             ,
	MinEventDuration,
	GapLength               ,
	NumSamplesForDerivativeEstimation,
	MedianValue             ,
	AvgValue                ,
	TruncateReal            ,
	TruncateHilbert         ,
	MAX_PARAM

}MotionXMLParams;

typedef struct
{
	void* parameter ;
}xmlParameterAddr ;


void parseDoc(char *docname , MotionXML_Params  *Params );
void XmlPrintStuctureParameters( MotionXML_Params *Params );


#endif /* BGXMLPARSER_H_ */



#ifdef __cplusplus
}
#endif
