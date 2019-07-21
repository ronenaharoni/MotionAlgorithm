#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <stddef.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <arpa/inet.h>
#include <sys/select.h>
#include <sys/socket.h>
#include <unistd.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <termios.h>
#include <netinet/in.h>
#include <time.h>
#include <sys/time.h>
#include <sys/mman.h> /* mmap() is defined in this header */
#include <libxml/parser.h>
#include <libxml/xmlmemory.h>
#include "BgXmlParser.h"
//#include "Utilities.h"

//#include <complex.h>

//_____________________________________________________________________________
//
// #defines
//_____________________________________________________________________________

//_____________________________________________________________________________
//
// typedefs
//_____________________________________________________________________________



//_____________________________________________________________________________
//
// static data
//_____________________________________________________________________________
#if 0
BgRadarParams bgParams ;

xmlParameterAddr parmChar[] =
{
		bgParams.IpAddr                   ,
		&bgParams.port                    ,
		&bgParams.debugFileFlag           ,
		&bgParams.Rstart_m                ,
		&bgParams.Rstop_m                 ,
		&bgParams.PII                     ,
		&bgParams.Tzero_ns                ,
		&bgParams.codeChannel             ,
		&bgParams.transmitGain            ,
		&bgParams.debugMode               ,
		&bgParams.Nscans_Resp             ,
		&bgParams.TxBitForMotion          ,
		&bgParams.RxAntForMotion          ,
		&bgParams.Rmin                    ,
		&bgParams.Acq_Nscans_Motion       ,
		&bgParams.EnableAcq               ,
		&bgParams.AcqDebugMode            ,
		&bgParams.saveFramesFlag          ,
		&bgParams.delta_R                 ,
		bgParams.DelayBetweenRecords      ,
};
#endif
//_____________________________________________________________________________
//
// Private function prototypes
//_____________________________________________________________________________
//_____________________________________________________________________________
//
// mrmSampleExit - close interfaces and exit sample app
//_____________________________________________________________________________

void XmlInitializeParamsAddressArray( xmlParameterAddr *paramAddresArray , MotionXML_Params *MotionXML_Params)
{
	int i = 0 ;

	(paramAddresArray[i++]).parameter = &MotionXML_Params->Nscans;
	(paramAddresArray[i++]).parameter = &MotionXML_Params->Nbins;
	(paramAddresArray[i++]).parameter = &MotionXML_Params->Fs;
	(paramAddresArray[i++]).parameter = &MotionXML_Params->RF_NumOfClasses;
	(paramAddresArray[i++]).parameter = &MotionXML_Params->SpectrogramRangeWidth;
	(paramAddresArray[i++]).parameter = &MotionXML_Params->SpectrogramWinLength;
	(paramAddresArray[i++]).parameter = &MotionXML_Params->SpectrogramN_Overlap;
	(paramAddresArray[i++]).parameter = &MotionXML_Params->SpectrogramNoiseFreqBins;
	(paramAddresArray[i++]).parameter = &MotionXML_Params->FminBin;
	(paramAddresArray[i++]).parameter = &MotionXML_Params->Fbins;
	(paramAddresArray[i++]).parameter = &MotionXML_Params->NoiseThresh;
	(paramAddresArray[i++]).parameter = &MotionXML_Params->TopMaxFreq;
	(paramAddresArray[i++]).parameter = &MotionXML_Params->MinFreqPlusMinus;
	(paramAddresArray[i++]).parameter = &MotionXML_Params->MinFreqReal;
	(paramAddresArray[i++]).parameter = &MotionXML_Params->MinEventDuration;
	(paramAddresArray[i++]).parameter = &MotionXML_Params->GapLength;
	(paramAddresArray[i++]).parameter = &MotionXML_Params->NumSamplesForDerivativeEstimation;
	(paramAddresArray[i++]).parameter = &MotionXML_Params->MedianValue;
	(paramAddresArray[i++]).parameter = &MotionXML_Params->AvgValue;
	(paramAddresArray[i++]).parameter = &MotionXML_Params->TruncateReal;
	(paramAddresArray[i++]).parameter = &MotionXML_Params->TruncateHilbert;


}








void XmlPrintStuctureParameters( MotionXML_Params *MotionXML_Params )
{
	xmlParameterAddr parmChar[TruncateHilbert] = { 0 } ;
	int index = 0 ;

	if( MotionXML_Params == NULL )
		return ;

	XmlInitializeParamsAddressArray( parmChar , MotionXML_Params );

	printf("print all parameters in structure \n");
	for( index = Nscans ; index <= TruncateHilbert ; index++ )
	{
		switch( index )
		{
		case Nscans             :
			///			printf( "%s\n" , (char*)parmChar[index].parameter );
			break;
		case Nbins :
		case Fs :
		case RF_NumOfClasses         :
		case SpectrogramRangeWidth   :
		case SpectrogramWinLength          :
			//			printf( "%d\n" , *(mrm_uint16_t*)(parmChar[index].parameter) );
			break;
		case SpectrogramN_Overlap     :
		case SpectrogramNoiseFreqBins                 :
		case FminBin            :
		case Fbins        :
		case NoiseThresh       :
		case TopMaxFreq          :
		case MinFreqPlusMinus   :
		case MinFreqReal   :
		case MinEventDuration                :
		case GapLength      :
		case NumSamplesForDerivativeEstimation    :
			//			printf( "%d\n" , *(mrm_uint8_t*)(parmChar[index].parameter) );
			break ;
		case MedianValue            :
		case AvgValue             :
		case TruncateReal             :
		case TruncateHilbert   :

		default :
			printf("out of range parameter \n");
			break;

		}
	}
}







void parseStory (xmlDocPtr doc, xmlNodePtr cur)
{
	xmlChar *key;
	cur = cur->xmlChildrenNode;
	while (cur != NULL)
	{
		if ((!xmlStrcmp(cur->name, (const xmlChar *)"keyword")))
		{
			key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
			printf("keyword: %s\n", key);
			xmlFree(key);
		}
		if ((!xmlStrcmp(cur->name, (const xmlChar *)"author")))
		{
			key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
			//			printf("author: %s\n", key);
			xmlFree(key);
		}
		cur = cur->next;
	}
	return;
}


void parse (xmlDocPtr doc, xmlNodePtr cur , MotionXML_Params *MotionXML_Params )
{
	xmlChar *key;
	int index = Nscans ;
	xmlParameterAddr parmChar[MAX_PARAM] = { 0 } ;

	if( cur == NULL )
		return ;

	XmlInitializeParamsAddressArray( parmChar , MotionXML_Params);

	cur = cur->xmlChildrenNode ;
	while(cur)
	{
		if( cur->type == XML_ELEMENT_NODE )
		{
			key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
			printf("%s: %s\n", cur->name , key );
			switch( index )
			{
			case Nbins :
				sscanf((const char*)key , "%d" , (int*)(parmChar[index].parameter) );
				break;
			case Nscans             :
				sscanf((const char*)key , "%d" , (int*)(parmChar[index].parameter) );
				break;
			case Fs :
				sscanf((const char*)key , "%d" , (int*)(parmChar[index].parameter) );
				break;
			case RF_NumOfClasses         :
				sscanf((const char*)key , "%d" , (int*)(parmChar[index].parameter) );
				break;
			case SpectrogramRangeWidth   :
				sscanf((const char*)key , "%d" , (int*)(parmChar[index].parameter) );
				break;
			case SpectrogramWinLength          :
				sscanf((const char*)key , "%d" , (int*)(parmChar[index].parameter) );
				break;
			case SpectrogramN_Overlap     :
				sscanf((const char*)key , "%d" , (int*)(parmChar[index].parameter) );
				break;
			case SpectrogramNoiseFreqBins                 :
				sscanf((const char*)key , "%d" , (int*)(parmChar[index].parameter) );
				break;
			case FminBin            :
				sscanf((const char*)key , "%d" , (int*)(parmChar[index].parameter) );
				break;
			case Fbins        :
				sscanf((const char*)key , "%d" , (int*)(parmChar[index].parameter) );
				break;
			case NoiseThresh       :
				sscanf((const char*)key , "%d" , (int*)(parmChar[index].parameter) );
				break;
			case TopMaxFreq          :
				sscanf((const char*)key , "%d" , (int*)(parmChar[index].parameter) );
				break;
			case MinFreqPlusMinus   :
				sscanf((const char*)key , "%d" , (int*)(parmChar[index].parameter) );
				break;
			case MinFreqReal   :
				sscanf((const char*)key , "%f" , (float*)(parmChar[index].parameter) );
				break;
			case MinEventDuration                :
				sscanf((const char*)key , "%d" , (int*)(parmChar[index].parameter) );
				break;
			case GapLength      :
				sscanf((const char*)key , "%d" , (int*)(parmChar[index].parameter) );
				break;
			case NumSamplesForDerivativeEstimation    :
				sscanf((const char*)key , "%d" , (int*)(parmChar[index].parameter) );
				break ;
			case MedianValue            :
				sscanf((const char*)key , "%d" , (int*)(parmChar[index].parameter) );
				break;
			case AvgValue             :
				sscanf((const char*)key , "%d" , (int*)(parmChar[index].parameter) );
				break;
			case TruncateReal             :
				sscanf((const char*)key , "%d" , (int*)(parmChar[index].parameter) );
				break;
			case TruncateHilbert   :
				sscanf((const char*)key , "%d" , (int*)(parmChar[index].parameter) );

				break;
//			case ROOM_COORDINATE:
//			{
//				float *pFloat = (float*) parmChar[index].parameter ;
//				StringReplaceChars( (char*)key , '[' , ' ' );
//				StringReplaceChars( (char*)key , ',' , ' ' );
//				StringReplaceChars( (char*)key , ']' , ' ' );
//
//				sscanf(((const char*)key) , "%f%f%f%f" , pFloat,pFloat+1,pFloat+2 , pFloat+3 );
//				break;
//			}
			default :
				printf("out of range parameter \n");
				break;
			}
			xmlFree(key);
			index++ ;


		}

#if 0
		//			if ((!xmlStrcmp(cur->name, (const xmlChar *)"port")))
		//			{
		if ( index == IP_ADDR )
		{
			printf("ip size %d %s" , (size_t)strlen( (const char*)key) , (char*)key  );
			memset( (char*)parmChar[index].parameter , 0 , 20 );
			memcpy( (char*)parmChar[index].parameter , (char*)key , (size_t)strlen( (const char*)key) );
		}
		if( index == IP_ADDR+1 )
		{
			sscanf((const char*)key , "%d" , (int*)(parmChar[index].parameter) );
			printf("%d --- %d \n" ,index ,  *(int*)(parmChar[index].parameter));

		}

		else if ( index == DELAY_BETWEEN_RECORDS )
		{
			float *pFloat = (float*) parmChar[index].parameter ;
			sscanf((const char*)key , "%f%f%f%f" , pFloat,pFloat+1,pFloat+2 , pFloat+3 );
		}
		else
		{
			sscanf((const char*)key , "%d" , (int*)(parmChar[index].parameter) );
			printf("%d --- %d \n" ,index ,  *(int*)(parmChar[index].parameter));

		}
#endif


		//printf("params: %s ---%d \n", key , Params->port );
		//			}


		cur = cur->next ;
	}

	//	XmlPrintStuctureParameters( Params );
#if 0
	printf("print all parameters in structure \n");
	for( index = IP_ADDR ; index <= DELAY_BETWEEN_RECORDS ; index++ )
	{
		switch( index )
		{
		case IP_ADDR:
			printf( "%s\n" , (char*)parmChar[index].parameter );
			break;
		case PORT :
		case NSCANS_RESP :
		case ACQ_NSCANS_MOTION :
		case ENABLE_ACQ :
			printf( "%d\n" , *(mrm_uint16_t*)(parmChar[index].parameter) );
			break;
		case DEBUG_FILE_FLAG  :
		case PII              :
		case TZERO_NS         :
		case CODE_CHANNEL     :
		case TRANSMIT_GAIN    :
		case DEBUG_MODE       :
		case TX_BIT_FOR_MOTION:
		case RX_ANT_FOR_MOTION:
		case RMIN             :
		case ACQ_DEBUG_MODE   :
		case SAVE_FRAMES_FLAG :
			printf( "%d\n" , *(mrm_uint8_t*)(parmChar[index].parameter) );
			break ;
		case RSTART_M         :
		case RSTOP_M          :
		case DELTA_R          :
			printf("%f\n", *(float*)(parmChar[index].parameter));
			break;
		case DELAY_BETWEEN_RECORDS :
			printf( "%f %f %f %f " , ((float*)(parmChar[index].parameter))[0]
																		   , ((float*)(parmChar[index].parameter))[1]
																												   , ((float*)(parmChar[index].parameter))[2]
																																						   , ((float*)(parmChar[index].parameter))[3]);
			break;

		default :
			printf("out of range parameter \n");
			break;

		}
	}
#endif
}

void parseDoc(char *docname ,MotionXML_Params *MotionXML_Params)
{
	xmlDocPtr doc;
	xmlNodePtr cur;
	doc = xmlParseFile(docname);
	if (doc == NULL )
	{
		fprintf(stderr,"Document not parsed successfully. \n");
		return;
	}
	cur = xmlDocGetRootElement(doc);
	if (cur == NULL)
	{
		fprintf(stderr,"empty document\n");
		xmlFreeDoc(doc);
		return;
	}
	if (xmlStrcmp(cur->name, (const xmlChar *) "ECHO_Conf"))
	{
		fprintf(stderr,"document of the wrong type, root node != story");
		xmlFreeDoc(doc);
		return;
	}
	cur = cur->xmlChildrenNode;

	while (cur != NULL)
	{
		if ((!xmlStrcmp(cur->name, (const xmlChar *)"Motion")))
		{
			parse (doc, cur , MotionXML_Params );
		}
		cur = cur->next;
	}
	xmlFreeDoc(doc);
	return;
}



