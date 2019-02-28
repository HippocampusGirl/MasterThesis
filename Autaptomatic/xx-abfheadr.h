//Filename:       xx-abfheadr.h
//Author:         Pamela Y. C. Fong
//Maintainer:     Qi Huang
//PI:             Henry A. Lester
//Organization:   California Institute of Technology
//
//                Reformulated and adapted from API copyrighted
//                by Axon Instruments, Inc.
//
//Created:        February 27, 2004
//Last Modified:  Mon Mar  1 13:16:04 PST 2004
//
//Instructions for constructing this file. 
//  Copy lines 54-974 from the AxonDev 1.8 API file, ABFHEADR.H,
//  commenting out the ASSERT statement on line 961.
//  Add include and namespace statements
//PYCF

#include <iostream>
#include <fstream>
#include <string>
#include <sys/stat.h>
using namespace std;

//
// Constants used in defining the ABF file header
//

#define ABF_ADCCOUNT           16    // number of ADC channels supported.
#define ABF_DACCOUNT           4     // number of DAC channels supported.
#define ABF_WAVEFORMCOUNT      2     // number of DAC channels which support waveforms.
#define ABF_EPOCHCOUNT         10    // number of waveform epochs supported.
#define ABF_BELLCOUNT          2     // Number of auditory signals supported.
#define ABF_ADCUNITLEN         8     // length of ADC units strings
#define ABF_ADCNAMELEN         10    // length of ADC channel name strings
#define ABF_DACUNITLEN         8     // length of DAC units strings
#define ABF_DACNAMELEN         10    // length of DAC channel name strings
#define ABF_VARPARAMLISTLEN    80    // length of conditioning string
#define ABF_USERLISTLEN        256   // length of the user list (V1.6)
#define ABF_USERLISTCOUNT      4     // number of independent user lists (V1.6)
#define ABF_OLDFILECOMMENTLEN  56    // length of file comment string (pre V1.6)
#define ABF_FILECOMMENTLEN     128   // length of file comment string (V1.6)

#define ABF_CREATORINFOLEN     16    // length of file creator info string
#define ABF_OLDDACFILENAMELEN  12    // old length of the DACFile name string
#define ABF_OLDDACFILEPATHLEN  60    // old length of the DACFile path string
#define ABF_DACFILEPATHLEN     84    // length of full path for DACFile
#define ABF_PATHLEN            256   // length of full path, used for DACFile and Protocol name.
#define ABF_ARITHMETICOPLEN    2     // length of the Arithmetic operator field
#define ABF_ARITHMETICUNITSLEN 8     // length of arithmetic units string
#define ABF_TAGCOMMENTLEN      56    // length of tag comment string
#define ABF_LONGDESCRIPTIONLEN 56    // length of long description entry
#define ABF_NOTENAMELEN        10    // length of the name component of a note
#define ABF_NOTEVALUELEN       8     // length of the value component of a note
#define ABF_NOTEUNITSLEN       8     // length of the units component of a note
#define ABF_BLOCKSIZE          512   // Size of block alignment in ABF files.
#define ABF_MACRONAMELEN       64    // Size of a Clampfit macro name.

#define ABF_CURRENTVERSION     1.80F           // Current file format version number
#define ABF_PREVIOUSVERSION    1.5F            // Previous file format version number (for old header size)
#define ABF_V16                1.6F            // Version number when the header size changed.
#define ABF_HEADERSIZE         6144            // Size of a Version 1.6 or later header
#define ABF_OLDHEADERSIZE      2048            // Size of a Version 1.5 or earlier header
#define ABF_NATIVESIGNATURE    0x20464241      // PC="ABF ", MAC=" FBA"
#define ABF_REVERSESIGNATURE   0x41424620      // PC=" FBA", MAC="ABF "

#define PCLAMP6_MAXSWEEPLENGTH         16384   // Maximum multiplexed sweep length supported by pCLAMP6 apps.
#define PCLAMP7_MAXSWEEPLEN_PERCHAN    1032258  // Maximum per channel sweep length supported by pCLAMP7 apps.

#define ABF_MAX_TRIAL_SAMPLES  0x7FFFFFFF    // Maximum length of acquisition supported (samples)
                                             // INT_MAX is used instead of UINT_MAX because of the signed 
                                             // values in the ABF header.

#define ABF_MAX_SWEEPS_PER_AVERAGE 65500     // The maximum number of sweeps that can be combined into a
                                             // cumulative average (nAverageAlgorithm=ABF_INFINITEAVERAGE).

#define ABF_STATS_REGIONS     8              // The number of independent statistics regions.
#define ABF_BASELINE_REGIONS  1              // The number of independent baseline regions.

#ifdef _MAC
   #define ABF_OLDPCLAMP        ABF_REVERSESIGNATURE
#else
   #define ABF_OLDPCLAMP        ABF_NATIVESIGNATURE
#endif

//
// Constant definitions for nFileType
//
#define ABF_ABFFILE          1
#define ABF_FETCHEX          2
#define ABF_CLAMPEX          3

//
// Constant definitions for nDataFormat
//
#define ABF_INTEGERDATA      0
#define ABF_FLOATDATA        1

//
// Constant definitions for nOperationMode
//
#define ABF_VARLENEVENTS     1
#define ABF_FIXLENEVENTS     2     // (ABF_FIXLENEVENTS == ABF_LOSSFREEOSC)
#define ABF_LOSSFREEOSC      2
#define ABF_GAPFREEFILE      3
#define ABF_HIGHSPEEDOSC     4
#define ABF_WAVEFORMFILE     5

//
// Constant definitions for nParamToVary
//
#define ABF_CONDITNUMPULSES         0
#define ABF_CONDITBASELINEDURATION  1
#define ABF_CONDITBASELINELEVEL     2
#define ABF_CONDITSTEPDURATION      3
#define ABF_CONDITSTEPLEVEL         4
#define ABF_CONDITPOSTTRAINDURATION 5
#define ABF_CONDITPOSTTRAINLEVEL    6
#define ABF_EPISODESTARTTOSTART     7
#define ABF_INACTIVEHOLDING         8
#define ABF_DIGITALHOLDING          9
#define ABF_PNNUMPULSES             10
#define ABF_PARALLELVALUE           11
#define ABF_EPOCHINITLEVEL          (ABF_PARALLELVALUE + ABF_EPOCHCOUNT)
#define ABF_EPOCHINITDURATION       (ABF_EPOCHINITLEVEL + ABF_EPOCHCOUNT)
#define ABF_EPOCHTRAINPERIOD        (ABF_EPOCHINITDURATION + ABF_EPOCHCOUNT)
#define ABF_EPOCHTRAINPULSEWIDTH    (ABF_EPOCHTRAINPERIOD + ABF_EPOCHCOUNT)
// Next value is (ABF_EPOCHINITDURATION + ABF_EPOCHCOUNT)

//
// Constants for nAveragingMode
//
#define ABF_NOAVERAGING       0
#define ABF_SAVEAVERAGEONLY   1
#define ABF_AVERAGESAVEALL    2

//
// Constants for nAverageAlgorithm
//
#define ABF_INFINITEAVERAGE   0
#define ABF_SLIDINGAVERAGE    1

//
// Constants for nEpochType
//
#define ABF_EPOCHDISABLED           0     // disabled epoch
#define ABF_EPOCHSTEPPED            1     // stepped waveform
#define ABF_EPOCHRAMPED             2     // ramp waveform
#define ABF_EPOCH_TYPE_RECTANGLE    3     // rectangular pulse train
#define ABF_EPOCH_TYPE_TRIANGLE     4     // triangular waveform
#define ABF_EPOCH_TYPE_COSINE       5     // cosinusoidal waveform
#define ABF_EPOCH_TYPE_RESISTANCE   6     // resistance waveform
#define ABF_EPOCH_TYPE_BIPHASIC     7     // biphasic pulse train

//
// Constants for epoch resistance
//
#define ABF_MIN_EPOCH_RESISTANCE_DURATION 8

//
// Constants for nWaveformSource
//
#define ABF_WAVEFORMDISABLED     0               // disabled waveform
#define ABF_EPOCHTABLEWAVEFORM   1
#define ABF_DACFILEWAVEFORM      2

//
// Constants for nInterEpisodeLevel & nDigitalInterEpisode
//
#define ABF_INTEREPI_USEHOLDING    0
#define ABF_INTEREPI_USELASTEPOCH  1

//
// Constants for nExperimentType
//
#define ABF_VOLTAGECLAMP         0
#define ABF_CURRENTCLAMP         1
#define ABF_SIMPLEACQUISITION    2

//
// Constants for nAutosampleEnable
//
#define ABF_AUTOSAMPLEDISABLED   0
#define ABF_AUTOSAMPLEAUTOMATIC  1
#define ABF_AUTOSAMPLEMANUAL     2

//
// Constants for nAutosampleInstrument
//
#define ABF_INST_UNKNOWN         0   // Unknown instrument (manual or user defined telegraph table).
#define ABF_INST_AXOPATCH1       1   // Axopatch-1 with CV-4-1/100
#define ABF_INST_AXOPATCH1_1     2   // Axopatch-1 with CV-4-0.1/100
#define ABF_INST_AXOPATCH1B      3   // Axopatch-1B(inv.) CV-4-1/100
#define ABF_INST_AXOPATCH1B_1    4   // Axopatch-1B(inv) CV-4-0.1/100
#define ABF_INST_AXOPATCH201     5   // Axopatch 200 with CV 201
#define ABF_INST_AXOPATCH202     6   // Axopatch 200 with CV 202
#define ABF_INST_GENECLAMP       7   // GeneClamp
#define ABF_INST_DAGAN3900       8   // Dagan 3900
#define ABF_INST_DAGAN3900A      9   // Dagan 3900A
#define ABF_INST_DAGANCA1_1      10  // Dagan CA-1  Im=0.1
#define ABF_INST_DAGANCA1        11  // Dagan CA-1  Im=1.0
#define ABF_INST_DAGANCA10       12  // Dagan CA-1  Im=10
#define ABF_INST_WARNER_OC725    13  // Warner OC-725
#define ABF_INST_WARNER_OC725C   14  // Warner OC-725
#define ABF_INST_AXOPATCH200B    15  // Axopatch 200B
#define ABF_INST_DAGANPCONE0_1   16  // Dagan PC-ONE  Im=0.1
#define ABF_INST_DAGANPCONE1     17  // Dagan PC-ONE  Im=1.0
#define ABF_INST_DAGANPCONE10    18  // Dagan PC-ONE  Im=10
#define ABF_INST_DAGANPCONE100   19  // Dagan PC-ONE  Im=100
#define ABF_INST_WARNER_BC525C   20  // Warner BC-525C
#define ABF_INST_WARNER_PC505    21  // Warner PC-505
#define ABF_INST_WARNER_PC501    22  // Warner PC-501
#define ABF_INST_DAGANCA1_05     23  // Dagan CA-1  Im=0.05
#define ABF_INST_MULTICLAMP700   24  // MultiClamp 700
#define ABF_INST_TURBO_TEC       25  // Turbo Tec
#define ABF_INST_OPUSXPRESS6000  26  // OpusXpress 6000A

//
// Constants for nManualInfoStrategy
//
#define ABF_ENV_DONOTWRITE      0
#define ABF_ENV_WRITEEACHTRIAL  1
#define ABF_ENV_PROMPTEACHTRIAL 2

//
// Constants for nTriggerSource
//
#define ABF_TRIGGERLINEINPUT           -5   // Start on line trigger (DD1320 only)
#define ABF_TRIGGERTAGINPUT            -4
#define ABF_TRIGGERFIRSTCHANNEL        -3
#define ABF_TRIGGEREXTERNAL            -2
#define ABF_TRIGGERSPACEBAR            -1
// >=0 = ADC channel to trigger off.

//
// Constants for nTrialTriggerSource
//
#define ABF_TRIALTRIGGER_SWSTARTONLY   -6   // Start on software message, end when protocol ends.
#define ABF_TRIALTRIGGER_SWSTARTSTOP   -5   // Start and end on software messages.
#define ABF_TRIALTRIGGER_LINEINPUT     -4   // Start on line trigger (DD1320 only)
#define ABF_TRIALTRIGGER_SPACEBAR      -3   // Start on spacebar press.
#define ABF_TRIALTRIGGER_EXTERNAL      -2   // Start on external trigger high
#define ABF_TRIALTRIGGER_NONE          -1   // Start immediately (default).
// >=0 = ADC channel to trigger off.    // Not implemented as yet...

//
// Constants for nTriggerPolarity.
//
#define ABF_TRIGGER_RISINGEDGE  0
#define ABF_TRIGGER_FALLINGEDGE 1

//
// Constants for nTriggerAction
//
#define ABF_TRIGGER_STARTEPISODE 0
#define ABF_TRIGGER_STARTRUN     1
#define ABF_TRIGGER_STARTTRIAL   2    // N.B. Discontinued in favor of nTrialTriggerSource

//
// Constants for nDrawingStrategy
//
#define ABF_DRAW_NONE            0
#define ABF_DRAW_REALTIME        1
#define ABF_DRAW_FULLSCREEN      2
#define ABF_DRAW_ENDOFRUN        3

//
// Constants for nTiledDisplay
//
#define ABF_DISPLAY_SUPERIMPOSED 0
#define ABF_DISPLAY_TILED        1

//
// Constants for nDataDisplayMode
//
#define ABF_DRAW_POINTS       0
#define ABF_DRAW_LINES        1

//
// Constants for nArithmeticExpression
//
#define ABF_SIMPLE_EXPRESSION    0
#define ABF_RATIO_EXPRESSION     1

//
// Constants for nLowpassFilterType & nHighpassFilterType
//
#define ABF_FILTER_NONE          0
#define ABF_FILTER_EXTERNAL      1
#define ABF_FILTER_SIMPLE_RC     2
#define ABF_FILTER_BESSEL        3
#define ABF_FILTER_BUTTERWORTH   4

//
// Constants for nPNPosition
//
#define ABF_PN_BEFORE_EPISODE    0
#define ABF_PN_AFTER_EPISODE     1

//
// Constants for nPNPolarity
//
#define ABF_PN_OPPOSITE_POLARITY -1
#define ABF_PN_SAME_POLARITY     1

//
// Constants for nAutopeakPolarity
//
#define ABF_PEAK_NEGATIVE       -1
#define ABF_PEAK_ABSOLUTE        0
#define ABF_PEAK_POSITIVE        1

//
// Constants for nAutopeakSearchMode
//
#define ABF_PEAK_SEARCH_SPECIFIED       -2
#define ABF_PEAK_SEARCH_ALL             -1
// nAutopeakSearchMode 0..9   = epoch in waveform 0's epoch table
// nAutopeakSearchMode 10..19 = epoch in waveform 1's epoch table

//
// Constants for nAutopeakBaseline
//
#define ABF_PEAK_BASELINE_SPECIFIED    -3
#define ABF_PEAK_BASELINE_NONE 	      -2
#define ABF_PEAK_BASELINE_FIRSTHOLDING -1
#define ABF_PEAK_BASELINE_LASTHOLDING  -4

//
// Constants for lAutopeakMeasurements
//
#define ABF_PEAK_MEASURE_PEAK                0x00000001
#define ABF_PEAK_MEASURE_PEAKTIME            0x00000002
#define ABF_PEAK_MEASURE_ANTIPEAK            0x00000004
#define ABF_PEAK_MEASURE_ANTIPEAKTIME        0x00000008
#define ABF_PEAK_MEASURE_MEAN                0x00000010
#define ABF_PEAK_MEASURE_STDDEV              0x00000020
#define ABF_PEAK_MEASURE_INTEGRAL            0x00000040
#define ABF_PEAK_MEASURE_MAXRISESLOPE        0x00000080
#define ABF_PEAK_MEASURE_MAXRISESLOPETIME    0x00000100
#define ABF_PEAK_MEASURE_MAXDECAYSLOPE       0x00000200
#define ABF_PEAK_MEASURE_MAXDECAYSLOPETIME   0x00000400
#define ABF_PEAK_MEASURE_RISETIME            0x00000800
#define ABF_PEAK_MEASURE_DECAYTIME           0x00001000
#define ABF_PEAK_MEASURE_HALFWIDTH           0x00002000
#define ABF_PEAK_MEASURE_BASELINE            0x00004000
#define ABF_PEAK_MEASURE_RISESLOPE           0x00008000
#define ABF_PEAK_MEASURE_DECAYSLOPE          0x00010000
#define ABF_PEAK_MEASURE_REGIONSLOPE         0x00020000
#define ABF_PEAK_MEASURE_ALL                 0x0002FFFF    // All of the above OR'd together.

//
// Constants for nStatsActiveChannels
//
#define ABF_PEAK_SEARCH_CHANNEL0          0x0001
#define ABF_PEAK_SEARCH_CHANNEL1          0x0002
#define ABF_PEAK_SEARCH_CHANNEL2          0x0004
#define ABF_PEAK_SEARCH_CHANNEL3          0x0008
#define ABF_PEAK_SEARCH_CHANNEL4          0x0010
#define ABF_PEAK_SEARCH_CHANNEL5          0x0020
#define ABF_PEAK_SEARCH_CHANNEL6          0x0040
#define ABF_PEAK_SEARCH_CHANNEL7          0x0080
#define ABF_PEAK_SEARCH_CHANNEL8          0x0100
#define ABF_PEAK_SEARCH_CHANNEL9          0x0200
#define ABF_PEAK_SEARCH_CHANNEL10         0x0400
#define ABF_PEAK_SEARCH_CHANNEL11         0x0800
#define ABF_PEAK_SEARCH_CHANNEL12         0x1000
#define ABF_PEAK_SEARCH_CHANNEL13         0x2000
#define ABF_PEAK_SEARCH_CHANNEL14         0x4000
#define ABF_PEAK_SEARCH_CHANNEL15         0x8000
#define ABF_PEAK_SEARCH_CHANNELSALL       0xFFFF      // All of the above OR'd together.

// Bit flag settings for nStatsSearchRegionFlags
//
#define ABF_PEAK_SEARCH_REGION0           0x01
#define ABF_PEAK_SEARCH_REGION1           0x02
#define ABF_PEAK_SEARCH_REGION2           0x04
#define ABF_PEAK_SEARCH_REGION3           0x08
#define ABF_PEAK_SEARCH_REGION4           0x10
#define ABF_PEAK_SEARCH_REGION5           0x20
#define ABF_PEAK_SEARCH_REGION6           0x40
#define ABF_PEAK_SEARCH_REGION7           0x80
#define ABF_PEAK_SEARCH_REGIONALL         0xFF        // All of the above OR'd together.

//
// Constants for lStatisticsMeasurements
//
#define ABF_STATISTICS_ABOVETHRESHOLD     0x00000001
#define ABF_STATISTICS_EVENTFREQUENCY     0x00000002
#define ABF_STATISTICS_MEANOPENTIME       0x00000004
#define ABF_STATISTICS_MEANCLOSEDTIME     0x00000008
#define ABF_STATISTICS_ALL                0x0000000F     // All the above OR'd together.

//
// Constants for nStatisticsSaveStrategy
//
#define ABF_STATISTICS_NOAUTOSAVE            0
#define ABF_STATISTICS_AUTOSAVE              1
#define ABF_STATISTICS_AUTOSAVE_AUTOCLEAR    2

//
// Constants for nStatisticsDisplayStrategy
//
#define ABF_STATISTICS_DISPLAY      0
#define ABF_STATISTICS_NODISPLAY    1

//
// Constants for nStatisticsClearStrategy
// determines whether to clear statistics after saving.
//
#define ABF_STATISTICS_NOCLEAR      0
#define ABF_STATISTICS_CLEAR        1

//
// Constants for nDACFileEpisodeNum
//
#define ABF_DACFILE_SKIPFIRSTSWEEP -1
#define ABF_DACFILE_USEALLSWEEPS    0
// >0 = The specific sweep number.

//
// Constants for nUndoPromptStrategy
//
#define ABF_UNDOPROMPT_ONABORT   0
#define ABF_UNDOPROMPT_ALWAYS    1

//
// Constants for nAutoAnalyseEnable
//
#define ABF_AUTOANALYSE_DISABLED   0
#define ABF_AUTOANALYSE_DEFAULT    1
#define ABF_AUTOANALYSE_RUNMACRO   2

//
// Constants for post nPostprocessLowpassFilterType
//
#define ABF_POSTPROCESS_FILTER_NONE          0
#define ABF_POSTPROCESS_FILTER_COMBINATION   1
#define ABF_POSTPROCESS_FILTER_BESSEL        2
#define ABF_POSTPROCESS_FILTER_BOXCAR        3
#define ABF_POSTPROCESS_FILTER_BUTTERWORTH   4
#define ABF_POSTPROCESS_FILTER_CHEBYSHEV     5
#define ABF_POSTPROCESS_FILTER_GAUSSIAN      6
#define ABF_POSTPROCESS_FILTER_RC            7
#define ABF_POSTPROCESS_FILTER_RC8           8
#define ABF_POSTPROCESS_FILTER_ADAPTIVE      9

//
// Miscellaneous constants
//
#define ABF_FILTERDISABLED  100000.0F     // Large frequency to disable lowpass filters
#define ABF_UNUSED_CHANNEL  -1            // Unused ADC and DAC channels.

//
// The output sampling sequence identifier for a seperate digital out channel.
//
#define ABF_DIGITAL_OUT_CHANNEL -1
#define ABF_PADDING_OUT_CHANNEL -2

//
// maximum values for various parameters (used by ABFH_CheckUserList).
//
#define ABF_CTPULSECOUNT_MAX           10000
#define ABF_CTBASELINEDURATION_MAX     100000.0F
#define ABF_CTSTEPDURATION_MAX         100000.0F
#define ABF_CTPOSTTRAINDURATION_MAX    100000.0F
#define ABF_SWEEPSTARTTOSTARTTIME_MAX  100000.0F 
#define ABF_PNPULSECOUNT_MAX           8
#define ABF_DIGITALVALUE_MAX           0xFF
#define ABF_EPOCHDIGITALVALUE_MAX      0x0F

//
// LTP Types - Reflects whether the header is used for LTP as baseline or induction.
//
#define ABF_LTP_TYPE_NONE              0
#define ABF_LTP_TYPE_BASELINE          1
#define ABF_LTP_TYPE_INDUCTION         2

//
// LTP Usage of DAC - Reflects whether the analog output will be used presynaptically or postsynaptically.
//
#define ABF_LTP_DAC_USAGE_NONE         0
#define ABF_LTP_DAC_USAGE_PRESYNAPTIC  1
#define ABF_LTP_DAC_USAGE_POSTSYNAPTIC 2

//
// Header Version Numbers
//
#define ABF_V166  1.66F
#define ABF_V167  1.67F
#define ABF_V168  1.68F
#define ABF_V169  1.69F
#define ABF_V170  1.70F
#define ABF_V171  1.71F
#define ABF_V172  1.72F
#define ABF_V173  1.73F
#define ABF_V174  1.74F
#define ABF_V175  1.75F
#define ABF_V176  1.76F
#define ABF_V177  1.77F
#define ABF_V178  1.78F
#define ABF_V179  1.79F
#define ABF_V180  1.80F

//
// pack structure on byte boundaries
//
#ifndef RC_INVOKED
#pragma pack(push, 1)
#endif

//
// Definition of the ABF header structure.
//

struct ABFFileHeader               // The total header length = 6144 bytes.
{
public:
   // GROUP #1 - File ID and size information. (40 bytes)
   long     lFileSignature;
   float    fFileVersionNumber;
   short    nOperationMode;
   long     lActualAcqLength;
   short    nNumPointsIgnored;
   long     lActualEpisodes;
   long     lFileStartDate;         // YYYYMMDD
   long     lFileStartTime;
   long     lStopwatchTime;
   float    fHeaderVersionNumber;
   short    nFileType;
   short    nMSBinFormat;

   // GROUP #2 - File Structure (78 bytes)
   long     lDataSectionPtr;
   long     lTagSectionPtr;
   long     lNumTagEntries;
   long     lScopeConfigPtr;
   long     lNumScopes;
   long     _lDACFilePtr;
   long     _lDACFileNumEpisodes;
   char     sUnused001[4];
   long     lDeltaArrayPtr;
   long     lNumDeltas;
   long     lVoiceTagPtr;
   long     lVoiceTagEntries;
   long     lUnused002;
   long     lSynchArrayPtr;
   long     lSynchArraySize;
   short    nDataFormat;
   short    nSimultaneousScan;
   long     lStatisticsConfigPtr;
   long     lAnnotationSectionPtr;
   long     lNumAnnotations;
   char     sUnused003[2];

   // GROUP #3 - Trial hierarchy information (82 bytes)
   /** 
   The number of input channels we acquired.
   Do not access directly - use CABFHeader::get_channel_count_acquired
   */
   short    channel_count_acquired;

   /** 
   The number of input channels we recorded.
   Do not access directly - use CABFHeader::get_channel_count_recorded
   */
   short    nADCNumChannels;
   float    fADCSampleInterval;
      /*{{
      The documentation says these two sample intervals are the interval between multiplexed samples, but not all digitisers work like that.
      Instead, these are the per-channel sample rate divided by the number of channels.
      If the user chose 100uS and has two channels, this value will be 50uS.
      }}*/
   float    fADCSecondSampleInterval;
      /*{{
      // The two sample intervals must be an integer multiple (or submultiple) of each other.
      if (fADCSampleInterval > fADCSecondSampleInterval)
         ASSERT(fmod(fADCSampleInterval, fADCSecondSampleInterval) == 0.0);
      if (fADCSecondSampleInterval, fADCSampleInterval)
         ASSERT(fmod(fADCSecondSampleInterval, fADCSampleInterval) == 0.0);
      }}*/
   float    fSynchTimeUnit;
   float    fSecondsPerRun;

   /**
   * The total number of samples per episode, for the recorded channels only.
   * This does not include channels which are acquired but not recorded.
   *
   * This is the number of samples per episode per channel, times the number of recorded channels.
   *
   * If you want the samples per episode for one channel, you must divide this by get_channel_count_recorded().
   */
   long     lNumSamplesPerEpisode;
   long     lPreTriggerSamples;
   long     lEpisodesPerRun;
   long     lRunsPerTrial;
   long     lNumberOfTrials;
   short    nAveragingMode;
   short    nUndoRunCount;
   short    nFirstEpisodeInRun;
   float    fTriggerThreshold;
   short    nTriggerSource;
   short    nTriggerAction;
   short    nTriggerPolarity;
   float    fScopeOutputInterval;
   float    fEpisodeStartToStart;
   float    fRunStartToStart;
   float    fTrialStartToStart;
   long     lAverageCount;
   long     lClockChange;
   short    nAutoTriggerStrategy;

   // GROUP #4 - Display Parameters (44 bytes)
   short    nDrawingStrategy;
   short    nTiledDisplay;
   short    nEraseStrategy;           // N.B. Discontinued. Use scope config entry instead.
   short    nDataDisplayMode;
   long     lDisplayAverageUpdate;
   short    nChannelStatsStrategy;
   long     lCalculationPeriod;       // N.B. Discontinued. Use fStatisticsPeriod.
   long     lSamplesPerTrace;
   long     lStartDisplayNum;
   long     lFinishDisplayNum;
   short    nMultiColor;
   short    nShowPNRawData;
   float    fStatisticsPeriod;
   long     lStatisticsMeasurements;
   short    nStatisticsSaveStrategy;

   // GROUP #5 - Hardware information (16 bytes)
   float    fADCRange;
   float    fDACRange;
   long     lADCResolution;
   long     lDACResolution;

   // GROUP #6 Environmental Information (118 bytes)
   short    nExperimentType;
   short    _nAutosampleEnable;
   short    _nAutosampleADCNum;
   short    _nAutosampleInstrument;
   float    _fAutosampleAdditGain;
   float    _fAutosampleFilter;
   float    _fAutosampleMembraneCap;
   short    nManualInfoStrategy;
   float    fCellID1;
   float    fCellID2;
   float    fCellID3;
   char     sCreatorInfo[ABF_CREATORINFOLEN];
   char     _sFileComment[ABF_OLDFILECOMMENTLEN];
   short    nFileStartMillisecs;    // Milliseconds portion of lFileStartTime
   short    nCommentsEnable;
   char     sUnused003a[8];

   // GROUP #7 - Multi-channel information (1044 bytes)
   short    nADCPtoLChannelMap[ABF_ADCCOUNT];
   short    nADCSamplingSeq[ABF_ADCCOUNT];
   char     sADCChannelName[ABF_ADCCOUNT][ABF_ADCNAMELEN];
   char     sADCUnits[ABF_ADCCOUNT][ABF_ADCUNITLEN];
   float    fADCProgrammableGain[ABF_ADCCOUNT];
   float    fADCDisplayAmplification[ABF_ADCCOUNT];
   float    fADCDisplayOffset[ABF_ADCCOUNT];       
   float    fInstrumentScaleFactor[ABF_ADCCOUNT];  
   float    fInstrumentOffset[ABF_ADCCOUNT];       
   float    fSignalGain[ABF_ADCCOUNT];
   float    fSignalOffset[ABF_ADCCOUNT];
   float    fSignalLowpassFilter[ABF_ADCCOUNT];
   float    fSignalHighpassFilter[ABF_ADCCOUNT];
   char     sDACChannelName[ABF_DACCOUNT][ABF_DACNAMELEN];
   char     sDACChannelUnits[ABF_DACCOUNT][ABF_DACUNITLEN];
   float    fDACScaleFactor[ABF_DACCOUNT];
   float    fDACHoldingLevel[ABF_DACCOUNT];
   short    nSignalType;
   char     sUnused004[10];

   // GROUP #8 - Synchronous timer outputs (14 bytes)
   short    nOUTEnable;
   short    nSampleNumberOUT1;
   short    nSampleNumberOUT2;
   short    nFirstEpisodeOUT;
   short    nLastEpisodeOUT;
   short    nPulseSamplesOUT1;
   short    nPulseSamplesOUT2;

   // GROUP #9 - Epoch Waveform and Pulses (184 bytes)
   short    nDigitalEnable;
   short    _nWaveformSource;
   short    nActiveDACChannel;
   short    _nInterEpisodeLevel;
   short    _nEpochType[ABF_EPOCHCOUNT];
   float    _fEpochInitLevel[ABF_EPOCHCOUNT];
   float    _fEpochLevelInc[ABF_EPOCHCOUNT];
   short    _nEpochInitDuration[ABF_EPOCHCOUNT];
   short    _nEpochDurationInc[ABF_EPOCHCOUNT];
   short    nDigitalHolding;
   short    nDigitalInterEpisode;
   short    nDigitalValue[ABF_EPOCHCOUNT];
   char     sUnavailable1608[4];    // was float fWaveformOffset;
   short    nDigitalDACChannel;
   char     sUnused005[6];

   // GROUP #10 - DAC Output File (98 bytes)
   float    _fDACFileScale;
   float    _fDACFileOffset;
   char     sUnused006[2];
   short    _nDACFileEpisodeNum;
   short    _nDACFileADCNum;
   char     _sDACFilePath[ABF_DACFILEPATHLEN];

   // GROUP #11 - Presweep (conditioning) pulse train (44 bytes)
   short    _nConditEnable;
   short    _nConditChannel;
   long     _lConditNumPulses;
   float    _fBaselineDuration;
   float    _fBaselineLevel;
   float    _fStepDuration;
   float    _fStepLevel;
   float    _fPostTrainPeriod;
   float    _fPostTrainLevel;
   char     sUnused007[12];

   // GROUP #12 - Variable parameter user list ( 82 bytes)
   short    _nParamToVary;
   char     _sParamValueList[ABF_VARPARAMLISTLEN];

   // GROUP #13 - Autopeak measurement (36 bytes)
   short    _nAutopeakEnable;
   short    _nAutopeakPolarity;
   short    _nAutopeakADCNum;
   short    _nAutopeakSearchMode;
   long     _lAutopeakStart;
   long     _lAutopeakEnd;
   short    _nAutopeakSmoothing;
   short    _nAutopeakBaseline;
   short    _nAutopeakAverage;
   char     sUnavailable1866[2];     // Was nAutopeakSaveStrategy, use nStatisticsSaveStrategy
   long     _lAutopeakBaselineStart;
   long     _lAutopeakBaselineEnd;
   long     _lAutopeakMeasurements;

   // GROUP #14 - Channel Arithmetic (52 bytes)
   short    nArithmeticEnable;
   float    fArithmeticUpperLimit;
   float    fArithmeticLowerLimit;
   short    nArithmeticADCNumA;
   short    nArithmeticADCNumB;
   float    fArithmeticK1;
   float    fArithmeticK2;
   float    fArithmeticK3;
   float    fArithmeticK4;
   char     sArithmeticOperator[ABF_ARITHMETICOPLEN];
   char     sArithmeticUnits[ABF_ARITHMETICUNITSLEN];
   float    fArithmeticK5;
   float    fArithmeticK6;
   short    nArithmeticExpression;
   char     sUnused008[2];

   // GROUP #15 - On-line subtraction (34 bytes)
   short    _nPNEnable;
   short    nPNPosition;
   short    _nPNPolarity;
   short    nPNNumPulses;
   short    _nPNADCNum;
   float    _fPNHoldingLevel;
   float    fPNSettlingTime;
   float    fPNInterpulse;
   char     sUnused009[12];

   // GROUP #16 - Miscellaneous variables (82 bytes)
   short    _nListEnable;
   
   short    nBellEnable[ABF_BELLCOUNT];
   short    nBellLocation[ABF_BELLCOUNT];
   short    nBellRepetitions[ABF_BELLCOUNT];
   
   short    nLevelHysteresis;
   long     lTimeHysteresis;
   short    nAllowExternalTags;
   
   char     nLowpassFilterType[ABF_ADCCOUNT];
   char     nHighpassFilterType[ABF_ADCCOUNT];
   short    nAverageAlgorithm;
   float    fAverageWeighting;
   short    nUndoPromptStrategy;
   short    nTrialTriggerSource;
   short    nStatisticsDisplayStrategy;
   short    nExternalTagType;
   long     lHeaderSize;
   double   dFileDuration;
   short    nStatisticsClearStrategy;
   // Size of v1.5 header = 2048

   // Extra parameters in v1.6
   // EXTENDED GROUP #2 - File Structure (26 bytes)
   long     lDACFilePtr[ABF_WAVEFORMCOUNT];
   long     lDACFileNumEpisodes[ABF_WAVEFORMCOUNT];
   char     sUnused010[10];
   
   // EXTENDED GROUP #7 - Multi-channel information (62 bytes)
   float    fDACCalibrationFactor[ABF_DACCOUNT];
   float    fDACCalibrationOffset[ABF_DACCOUNT];
   char     sUnused011[30];

   // GROUP #17 - Trains parameters (160 bytes)
   long     lEpochPulsePeriod[ABF_WAVEFORMCOUNT][ABF_EPOCHCOUNT];
   long     lEpochPulseWidth [ABF_WAVEFORMCOUNT][ABF_EPOCHCOUNT];

   // EXTENDED GROUP #9 - Epoch Waveform and Pulses ( 412 bytes)
   short    nWaveformEnable[ABF_WAVEFORMCOUNT];
   short    nWaveformSource[ABF_WAVEFORMCOUNT];
   short    nInterEpisodeLevel[ABF_WAVEFORMCOUNT];
   short    nEpochType[ABF_WAVEFORMCOUNT][ABF_EPOCHCOUNT];
   float    fEpochInitLevel[ABF_WAVEFORMCOUNT][ABF_EPOCHCOUNT];
   float    fEpochLevelInc[ABF_WAVEFORMCOUNT][ABF_EPOCHCOUNT];
   long     lEpochInitDuration[ABF_WAVEFORMCOUNT][ABF_EPOCHCOUNT];
   long     lEpochDurationInc[ABF_WAVEFORMCOUNT][ABF_EPOCHCOUNT];
   short    nDigitalTrainValue[ABF_EPOCHCOUNT];                         // 2 * 10 = 20 bytes
   short    nDigitalTrainActiveLogic;                                   // 2 bytes
   char     sUnused012[18];

   // EXTENDED GROUP #10 - DAC Output File (552 bytes)
   float    fDACFileScale[ABF_WAVEFORMCOUNT];
   float    fDACFileOffset[ABF_WAVEFORMCOUNT];
   long     lDACFileEpisodeNum[ABF_WAVEFORMCOUNT];
   short    nDACFileADCNum[ABF_WAVEFORMCOUNT];
   char     sDACFilePath[ABF_WAVEFORMCOUNT][ABF_PATHLEN];
   char     sUnused013[12];

   // EXTENDED GROUP #11 - Presweep (conditioning) pulse train (100 bytes)
   short    nConditEnable[ABF_WAVEFORMCOUNT];
   long     lConditNumPulses[ABF_WAVEFORMCOUNT];
   float    fBaselineDuration[ABF_WAVEFORMCOUNT];
   float    fBaselineLevel[ABF_WAVEFORMCOUNT];
   float    fStepDuration[ABF_WAVEFORMCOUNT];
   float    fStepLevel[ABF_WAVEFORMCOUNT];
   float    fPostTrainPeriod[ABF_WAVEFORMCOUNT];
   float    fPostTrainLevel[ABF_WAVEFORMCOUNT];
   char     sUnused014[40];

   // EXTENDED GROUP #12 - Variable parameter user list (1096 bytes)
   short    nULEnable[ABF_USERLISTCOUNT];
   short    nULParamToVary[ABF_USERLISTCOUNT];
   char     sULParamValueList[ABF_USERLISTCOUNT][ABF_USERLISTLEN];
   short    nULRepeat[ABF_USERLISTCOUNT];
   char     sUnused015[48];

   // EXTENDED GROUP #15 - On-line subtraction (56 bytes)
   short    nPNEnable[ABF_WAVEFORMCOUNT];
   short    nPNPolarity[ABF_WAVEFORMCOUNT];
   short    nPNADCNum[ABF_WAVEFORMCOUNT];
   float    fPNHoldingLevel[ABF_WAVEFORMCOUNT];
   char     sUnused016[36];

   // EXTENDED GROUP #6 Environmental Information  (898 bytes)
   short    nTelegraphEnable[ABF_ADCCOUNT];
   short    nTelegraphInstrument[ABF_ADCCOUNT];
   float    fTelegraphAdditGain[ABF_ADCCOUNT];
   float    fTelegraphFilter[ABF_ADCCOUNT];
   float    fTelegraphMembraneCap[ABF_ADCCOUNT];
   short    nTelegraphMode[ABF_ADCCOUNT];
   short    nTelegraphDACScaleFactorEnable[ABF_DACCOUNT];
   char     sUnused016a[24];

   short    nAutoAnalyseEnable;
   char     sAutoAnalysisMacroName[ABF_MACRONAMELEN];
   char     sProtocolPath[ABF_PATHLEN];

   char     sFileComment[ABF_FILECOMMENTLEN];
   char     sUnused017[128];

   // EXTENDED GROUP #13 - Statistics measurements (388 bytes)
   short    nStatsEnable;
   unsigned short nStatsActiveChannels;             // Active stats channel bit flag
   unsigned short nStatsSearchRegionFlags;          // Active stats region bit flag
   short    nStatsSelectedRegion;
   short    _nStatsSearchMode;
   short    nStatsSmoothing;
   short    nStatsSmoothingEnable;
   short    nStatsBaseline;
   long     lStatsBaselineStart;
   long     lStatsBaselineEnd;
   long     lStatsMeasurements[ABF_STATS_REGIONS];  // Measurement bit flag for each region
   long     lStatsStart[ABF_STATS_REGIONS];
   long     lStatsEnd[ABF_STATS_REGIONS];
   short    nRiseBottomPercentile[ABF_STATS_REGIONS];
   short    nRiseTopPercentile[ABF_STATS_REGIONS];
   short    nDecayBottomPercentile[ABF_STATS_REGIONS];
   short    nDecayTopPercentile[ABF_STATS_REGIONS];
   short    nStatsChannelPolarity[ABF_ADCCOUNT];
   short    nStatsSearchMode[ABF_STATS_REGIONS];    // Stats mode per region: mode is cursor region, epoch etc 
   char     sUnused018[156];

   // GROUP #18 - Application version data (16 bytes)
   short    nMajorVersion;
   short    nMinorVersion;
   short    nBugfixVersion;
   short    nBuildVersion;
   char     sUnused019[8];

   // GROUP #19 - LTP protocol (14 bytes)
   short    nLTPType;
   short    nLTPUsageOfDAC[ABF_WAVEFORMCOUNT];
   short    nLTPPresynapticPulses[ABF_WAVEFORMCOUNT];
   char     sUnused020[4];

   // GROUP #20 - Digidata 132x Trigger out flag. (8 bytes)
   short    nDD132xTriggerOut;
   char     sUnused021[6];

   // GROUP #21 - Epoch resistance (40 bytes)
   char     sEpochResistanceSignalName[ABF_WAVEFORMCOUNT][ABF_ADCNAMELEN];
   short    nEpochResistanceState[ABF_WAVEFORMCOUNT];
   char     sUnused022[16];
   
   // GROUP #22 - Alternating episodic mode (58 bytes)
   short    nAlternateDACOutputState;
   short    nAlternateDigitalValue[ABF_EPOCHCOUNT];
   short    nAlternateDigitalTrainValue[ABF_EPOCHCOUNT];
   short    nAlternateDigitalOutputState;
   char     sUnused023[14];

   // GROUP #23 - Post-processing actions (210 bytes)
   float    fPostProcessLowpassFilter[ABF_ADCCOUNT];
   char     nPostProcessLowpassFilterType[ABF_ADCCOUNT];


   // 6014 header bytes allocated + 130 header bytes not allocated
   char     sUnused2048[130];

   ABFFileHeader();
};   // Size = 6144
// This structure is persisted, so the size MUST NOT CHANGE
//STATIC_ASSERT(sizeof(ABFFileHeader) == 6144);

inline ABFFileHeader::ABFFileHeader()
{
   // Set everything to 0.
   memset( this, 0, sizeof(ABFFileHeader) );
   
   // Set critical parameters so we can determine the version.
   lFileSignature       = ABF_NATIVESIGNATURE;
   fFileVersionNumber   = ABF_CURRENTVERSION;
   fHeaderVersionNumber = ABF_CURRENTVERSION;
   lHeaderSize          = ABF_HEADERSIZE;
}

void ABFH_GetADCtoUUFactors(const ABFFileHeader *pFH, int nChannel, 
                                    float *pfADCToUUFactor, float *pfADCToUUShift );


