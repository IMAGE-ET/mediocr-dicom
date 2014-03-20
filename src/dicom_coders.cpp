#include <mediocr-dicom/dicom_coders.hpp>

#include <fstream>
#include <stdexcept>

#include <dcmtk/dcmdata/dcuid.h>
#include <dcmtk/dcmdata/dcmetinf.h>
#include <dcmtk/dcmdata/dcfilefo.h>
#include <dcmtk/dcmdata/dcdict.h>
#include <dcmtk/dcmdata/dcdeftag.h>
#include <dcmtk/dcmjpeg/djdecode.h>
#include <dcmtk/dcmjpeg/djencode.h>
#include <dcmtk/dcmjpls/djencode.h>
#include <dcmtk/dcmjpls/djdecode.h>
#include <dcmtk/dcmdata/dcrleerg.h>
#include <dcmtk/dcmdata/dcrledrg.h>
#include <dcmtk/ofstd/ofconapp.h>

namespace mediocr {

bool file_has_dicom_magic(std::string filepath) {
    // DICOM files open with a preamble of 128 bytes, then
    // the magic text "DICM".
    std::ifstream fp(filepath);
    if(!fp) {
        throw std::runtime_error("Failed to open file: " + filepath);
    }
    char prefix[5];
    fp.seekg(128);
    fp.read(prefix, 4);
    prefix[4] = 0;
    if(!fp) {
        // less than 132 bytes in file
        return false;
    }
    return std::string(prefix) == "DICM";
}

namespace dicom {

Coders::Coders()
{
	djencoder_register();
	djdecoder_register();

	DJLSEncoderRegistration::registerCodecs();
	DJLSDecoderRegistration::registerCodecs();

	DcmRLEEncoderRegistration::registerCodecs();
	DcmRLEDecoderRegistration::registerCodecs();
}

Coders::~Coders(){
	djdecoder_cleanup();
	djencoder_cleanup();

	DJLSEncoderRegistration::cleanup();
	DJLSDecoderRegistration::cleanup();

	DcmRLEEncoderRegistration::cleanup();
	DcmRLEDecoderRegistration::cleanup();
}

void Coders::djencoder_register() {
	E_CompressionColorSpaceConversion opt_compCSconversion = ECC_lossyYCbCr;
	E_UIDCreation opt_uidcreation = EUC_never;     // UID creation is controled by HeaderProcessed

	OFBool opt_huffmanOptimize = OFTrue;
	OFCmdUnsignedInt opt_smoothing = 0;
	int opt_compressedBits = 0;                    // 0=auto, 8/12/16=force
	OFCmdUnsignedInt opt_fragmentSize = 0;         // 0=unlimited
	OFBool opt_createOffsetTable = OFTrue;
	E_SubSampling opt_sampleFactors = ESS_444;
	OFBool opt_useYBR422 = OFFalse;
	OFBool opt_secondarycapture = OFFalse;
	int opt_windowType = 0;                        // 0=no windowing, 1=Wi, 2=Wl, 3=Wm, 4=Wh, 5=Ww, 6=Wn, 7=Wr
	OFCmdUnsignedInt opt_windowParameter = 0;
	OFCmdFloat opt_windowCenter=0.0;
	OFCmdFloat opt_windowWidth=0.0;
	OFCmdUnsignedInt opt_roiLeft = 0;
	OFCmdUnsignedInt opt_roiTop = 0;
	OFCmdUnsignedInt opt_roiWidth = 0;
	OFCmdUnsignedInt opt_roiHeight = 0;
	OFBool opt_usePixelValues = OFTrue;
	OFBool opt_useModalityRescale = OFFalse;
	OFBool opt_acceptWrongPaletteTags = OFFalse;
	OFBool opt_acrNemaCompatibility = OFFalse;
	OFBool opt_trueLossless = OFTrue;  // Newly added in dcmtk 3.5.4. Since this option defaults to false we need to include all options...

	DJEncoderRegistration::registerCodecs (
		opt_compCSconversion,
		opt_uidcreation,
		opt_huffmanOptimize,
		opt_smoothing,
		opt_compressedBits,
		opt_fragmentSize,
		opt_createOffsetTable,
		opt_sampleFactors,
		opt_useYBR422,
		opt_secondarycapture,
		opt_windowType,
		opt_windowParameter,
		opt_windowCenter,
		opt_windowWidth,
		opt_roiLeft,
		opt_roiTop,
		opt_roiWidth,
		opt_roiHeight,
		opt_usePixelValues,
		opt_useModalityRescale,
		opt_acceptWrongPaletteTags,
		opt_acrNemaCompatibility,
		opt_trueLossless
	);
}

void Coders::djencoder_cleanup(){
	DJEncoderRegistration::cleanup();
}

void Coders::djdecoder_register(){
	DJDecoderRegistration::registerCodecs();
}

void Coders::djdecoder_cleanup(){
	DJDecoderRegistration::cleanup();
}

}}
