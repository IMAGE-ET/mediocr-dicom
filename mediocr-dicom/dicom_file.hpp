#pragma once

#include <vector>
#include <stdexcept>
#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdint.h>

#include <dcmtk/dcmdata/dcuid.h>
#include <dcmtk/dcmdata/dcfilefo.h>
#include <dcmtk/dcmdata/dcdict.h>
#include <dcmtk/dcmdata/dcdeftag.h>

namespace mediocr {

namespace dicom {

template <typename Parameters>
struct transfer_syntax {
	E_TransferSyntax syntax;
	Parameters parameters;
};

//! Helper function to deduct Parameters' type
template <typename Parameters>
transfer_syntax<Parameters> create_transfer_syntax(decltype(transfer_syntax<Parameters>::syntax) const& s, Parameters p) {
	return transfer_syntax<Parameters>{s, p};
}

}

struct dicom_key_error : public std::runtime_error {
	DcmTagKey key;

	static std::string string_form(std::string prepend, DcmTagKey k) {
		std::stringstream ss;
		ss << prepend << k;
		return ss.str();
	}

	virtual ~dicom_key_error() throw () {}

	dicom_key_error(std::string s, DcmTagKey k)
	: runtime_error(string_form(s + ": ", k)), key(k)
	{}
};

/** A wrapper for a DCMTK-DcmFileFormat object */
struct dicom_file {
	dicom_file() {}

	dicom_file(std::string filename){
		if(header.loadFile(filename.c_str()).bad()){
			throw std::runtime_error("DCMTK failed to open (" + filename + ")");
		}
	}
	
	dicom_file(DcmFileFormat dicom)
	: header(dicom)
	{
		if(header.error().bad()){
			throw std::runtime_error("Tried to initialize dicom_file with a bad dicom");
		}
	}

	DcmDataset& get_dataset(){ return *header.getDataset(); }

	uint16_t get_bytes_per_pixel() const {
		return std::ceil(get_bits_allocated()/float((sizeof(char)*8)));
	}

	void loadAllDataIntoMemory() const {
		header.loadAllDataIntoMemory();
	}

	struct TagIterator {
		/* begin */
		TagIterator(dicom_file &file) : dataset(0), object(0) {
			dataset = file.header.getDataset();
			object = dataset->nextInContainer(object);
		}

		/* end */
		TagIterator() : dataset(0), object(0) {}

		bool operator==(const TagIterator &it) const {
			return dataset == it.dataset;
		}

		bool operator!=(const TagIterator &it) const {
			return dataset != it.dataset;
		}

		void operator++() {
			if(dataset) {
				// nextInContainer(NULL) returns the first object
				// nextInContainer(lastObject) returns NULL
				object = dataset->nextInContainer(object);
				if(object == NULL) {
					dataset = NULL;
				}
			}
		}

		void operator++(int) {
			++(*this);
		}

		DcmElement &operator*() const {
			return getElement();
		}

		DcmElement *operator->() const {
			return getElementPtr();
		}

		DcmElement &getElement() const {
			return dynamic_cast<DcmElement&>(*object);
		}

		DcmElement *getElementPtr() const {
			return dynamic_cast<DcmElement*>(object);
		}

		DcmTag getTag() const {
			return getElement().getTag();
		}

		std::string getTagName() const {
			return getTag().getTagName();
		}

		std::string getVRName() const {
			return getTag().getVRName();
		}

		bool isUnknownVR() const {
			return getTag().isUnknownVR();
		}

		uint32_t getValueLength() const {
			return getElement().getLength();
		}

		std::string getValueString() const {
			OFString str;
			if(getElement().getOFStringArray(str).bad()) {
				throw std::runtime_error("Failed to get value as a string: tag " + getTagName() + ", VR " + getVRName());
			}
			return std::string(str.data(), str.size());
		}

	private:
		DcmDataset *dataset;
		DcmObject *object;
	};

	TagIterator begin() {
		return TagIterator(*this);
	}

	TagIterator end() {
		return TagIterator();
	}

	void dumpTags(std::ostream &out) {
		for(auto it = begin(); it != end(); ++it) {
			out << it.getTagName();
			if(!it->isLeaf() || it.getVRName() == "OB" || it.getVRName() == "OW") {
				out << " = complex" << std::endl;
			} else {
				out << " = [" << it.getValueString() << "]" << std::endl;
			}
		}
	}

#define DCM_TAG_UINT16(func, tag) \
	uint16_t get_ ## func () const { \
		return get_uint16(tag); \
	} \
	void set_ ## func (uint16_t v) { \
		set_uint16(tag, v); \
	}
#define DCM_TAG_INTSTRING(func, tag) \
	int32_t get_ ## func () const { \
		int32_t res; \
		std::stringstream ss; \
		ss << get_string(tag); \
		ss >> res; \
		if(!ss) { \
			throw std::runtime_error("Invalid non-number string given in " # tag); \
		} \
		return res; \
	} \
	void set_ ## func (int32_t v) { \
		set_string(tag, std::to_string(v)); \
	}
#define DCM_TAG_STRING(func, tag) \
	std::string get_ ## func () const { \
		return get_string(tag); \
	} \
	void set_ ## func (std::string v) { \
		set_string(tag, v); \
	}

#define DCM_TAG_DOUBLE_GET(func, tag)\
	double get_ ## func () const { \
		return get_float64(tag); \
	}

#define DCM_TAG_DOUBLE(func, tag) \
	DCM_TAG_DOUBLE_GET(func, tag) \
	void set_ ## func (double v) { \
		set_float64(tag, v); \
	}
#define DCM_TAG_DOUBLE_ARRAY(func, tag, i) \
	double get_ ## func () const { \
		return get_float64(tag, i); \
	} \
	void set_ ## func(double v) { \
		set_float64(tag, v, i); \
	}

	DCM_TAG_UINT16(number_of_rows, DCM_Rows);
	DCM_TAG_UINT16(number_of_columns, DCM_Columns);
	DCM_TAG_UINT16(bits_allocated, DCM_BitsAllocated);
	DCM_TAG_UINT16(bits_stored, DCM_BitsStored);
	DCM_TAG_UINT16(pixel_representation, DCM_PixelRepresentation);
	DCM_TAG_UINT16(samples_per_pixel, DCM_SamplesPerPixel);
	DCM_TAG_UINT16(high_bit, DCM_HighBit);

	DCM_TAG_INTSTRING(number_of_planes, DCM_NumberOfFrames);
	DCM_TAG_INTSTRING(relative_xray_exposure, DCM_RelativeXRayExposure);
	DCM_TAG_INTSTRING(exposure, DCM_Exposure);
	DCM_TAG_INTSTRING(exposure_in_uas, DCM_ExposureInuAs);

	DCM_TAG_STRING(voi_lut_function, DCM_VOILUTFunction);
	DCM_TAG_STRING(study_uid, DCM_StudyInstanceUID);
	DCM_TAG_STRING(series_uid, DCM_SeriesInstanceUID);
	DCM_TAG_STRING(modality, DCM_Modality);
	DCM_TAG_STRING(patient_orientation, DCM_PatientOrientation);
	DCM_TAG_STRING(protocol_name, DCM_ProtocolName);
	DCM_TAG_STRING(manufacturer, DCM_Manufacturer);
	DCM_TAG_STRING(presentation_intent_type, DCM_PresentationIntentType);
	DCM_TAG_STRING(patient_name, DCM_PatientName);
	DCM_TAG_STRING(study_date, DCM_StudyDate);
	DCM_TAG_STRING(patient_birth_date, DCM_PatientBirthDate);
	DCM_TAG_STRING(patient_age, DCM_PatientAge);
	DCM_TAG_STRING(patient_sex, DCM_PatientSex);
	DCM_TAG_STRING(view_name, DCM_ViewName);
	DCM_TAG_STRING(view_position, DCM_ViewPosition);
	DCM_TAG_STRING(image_laterality, DCM_ImageLaterality);
	DCM_TAG_STRING(manufacturer_model_name, DCM_ManufacturerModelName);
	DCM_TAG_STRING(detector_type, DCM_DetectorType);
	DCM_TAG_STRING(anode_target_material, DCM_AnodeTargetMaterial);
	DCM_TAG_STRING(filter_material, DCM_FilterMaterial);
	DCM_TAG_STRING(photometric_interpretation, DCM_PhotometricInterpretation);

	DCM_TAG_DOUBLE_ARRAY(pixel_row_spacing, DCM_PixelSpacing, 0);
	DCM_TAG_DOUBLE_ARRAY(pixel_column_spacing, DCM_PixelSpacing, 1);
	DCM_TAG_DOUBLE_ARRAY(imager_pixel_row_spacing, DCM_ImagerPixelSpacing, 0);
	DCM_TAG_DOUBLE_ARRAY(imager_pixel_column_spacing, DCM_ImagerPixelSpacing, 1);
	DCM_TAG_DOUBLE_ARRAY(nominal_scanned_pixel_row_spacing, DCM_NominalScannedPixelSpacing, 0);
	DCM_TAG_DOUBLE_ARRAY(nominal_scanned_pixel_column_spacing, DCM_NominalScannedPixelSpacing, 1);

	DCM_TAG_DOUBLE(magnetic_field_strength, DCM_MagneticFieldStrength);
	DCM_TAG_DOUBLE(slice_thickness, DCM_SliceThickness);
	DCM_TAG_DOUBLE(slice_location, DCM_SliceLocation);
	DCM_TAG_DOUBLE_GET(window_width, DCM_WindowWidth);
	DCM_TAG_DOUBLE_GET(window_center, DCM_WindowCenter);
	DCM_TAG_DOUBLE_GET(rescale_intercept, DCM_RescaleIntercept);
	DCM_TAG_DOUBLE_GET(rescale_slope, DCM_RescaleSlope);
	DCM_TAG_DOUBLE_GET(spacing_between_slices, DCM_SpacingBetweenSlices);
	DCM_TAG_DOUBLE_GET(body_part_thickness, DCM_BodyPartThickness);
	DCM_TAG_DOUBLE_GET(kvp, DCM_KVP);
	DCM_TAG_DOUBLE_GET(compression_force, DCM_CompressionForce);

	std::vector<double> get_image_orientation_patient() const {
		std::vector<double> x;
		const unsigned int dimensions = 6;
		x.reserve(dimensions);
		
		for(unsigned int i = 0; i < dimensions; ++i){
			x.push_back(get_float64(DCM_ImageOrientationPatient, i));
		}
		
		return x;
	}
	
	std::vector<double> get_image_position_patient() const {
		std::vector<double> x;
		const unsigned int dimensions = 3;
		x.reserve(dimensions);
		
		for(unsigned int i = 0; i < dimensions; ++i){
			x.push_back(get_float64(DCM_ImagePositionPatient, i));
		}
		
		return x;
	}
	
	
	DcmFileFormat const& get_original_header() const {
		return header;
	}
	
	void decompress() {
		// Decompress if image is JPEG-encapsulated
		E_TransferSyntax Xfer = header.getDataset()->getOriginalXfer();
		OFBool hasRep = header.getDataset()->hasRepresentation (Xfer , NULL);

		/* unknown representation or JPEG-encapsulated in any form */
		if ( (!hasRep) || (Xfer >= 4)) {
			OFCondition condDecom = header.getDataset()->chooseRepresentation (EXS_LittleEndianExplicit, NULL);   //decompress
			header.getDataset()->removeAllButCurrentRepresentations();

			if (condDecom.bad()) {
				throw std::runtime_error("Error while decompressing JPEG DICOM image");
			}
		}
	}
	
	//! Get pixel data. Call decompress() first, or this function will return a nullptr.
	template <typename T>
	T* get_pixel_data_as() const;

	template <typename T>
	void set_pixel_data(T const& data){
		set_pixel_data(std::begin(data), std::end(data));
	}

	template <typename T>
	void set_pixel_data(std::vector<T> const& data) {
		set_pixel_data(&data.front(), &data.back() + 1);
	}

	//! We need the following setters because they are DSs and we can't set DS with setFloat64, only GET them with getFloat64
	void set_pixel_spacing(float const row_spacing, float const column_spacing){
		set_decimal_string(DCM_PixelSpacing, {row_spacing, column_spacing});
	}

	void set_imager_pixel_spacing(float const row_spacing, float const column_spacing){
		set_decimal_string(DCM_ImagerPixelSpacing, {row_spacing, column_spacing});
	}

	void set_rescale_slope(float const x){
		set_decimal_string(DCM_RescaleSlope, {x});
	}

	void set_rescale_intercept(float const x){
		set_decimal_string(DCM_RescaleIntercept, {x});
	}

	void set_window_width(float const x){
		set_decimal_string(DCM_WindowWidth, {x});
	}

	void set_window_center(float const x){
		set_decimal_string(DCM_WindowCenter, {x});
	}

	void set_spacing_between_slices(float const x){
		set_decimal_string(DCM_SpacingBetweenSlices, {x});
	}

	void set_body_part_thickness(float const x){
		set_decimal_string(DCM_BodyPartThickness, {x});
	}

	void set_kvp(float const x) {
		set_decimal_string(DCM_KVP, {x});
	}

	void set_image_position_patient(const std::vector<float> values) {
		set_decimal_string(DCM_ImagePositionPatient, values);
	}

	void set_image_orientation_patient(const std::vector<float> values) {
		set_decimal_string(DCM_ImageOrientationPatient, values);
	}

	void set_compression_force(float const x){
		set_decimal_string(DCM_CompressionForce, {x});
	}

	//! Saves a file in lossless JPEG encoding, best suited for mammograms (R. Visser, Lossless Compression of Digital Mammograms, IWDM 2006)
	void write_to_file(std::string filename);

private:
	void write_to_file(std::string filename, E_TransferSyntax transfer_syntax, DcmRepresentationParameter* p);

public:
	//! Saves a file in whatever syntax and parameters you supply. Check the DCMTK manual for more information.
	template <typename Parameters>
	void write_to_file(std::string filename, dicom::transfer_syntax<Parameters> t){
		write_to_file(filename, t.syntax, &t.parameters);
	}

	void remove_all(DcmTagKey key){
		remove(key, true);
	}

	void remove(DcmTagKey key, bool remove_all){
		if(header.getDataset()->findAndDeleteElement(key, remove_all).bad()){
			throw mediocr::dicom_key_error("Failed to remove some key", key);
		}
	}

	//! Anonymizes the DICOM by removing several tags as per 'DICOM Supplement 55: Attribute level confidentiality (including de-identification)'
	void anonymize() {
		static const DcmTagKey remove_these [] = {
			DCM_InstanceCreatorUID,
			DCM_SOPInstanceUID,
			DCM_AccessionNumber,
			DCM_InstitutionName,
			DCM_InstitutionAddress,
			DCM_ReferringPhysicianName,
			DCM_ReferringPhysicianAddress,
			DCM_ReferringPhysicianTelephoneNumbers,
			DCM_StationName,
			DCM_StudyDescription,
			DCM_SeriesDescription,
			DCM_InstitutionalDepartmentName,
			DCM_PhysiciansOfRecord,
			DCM_PerformingPhysicianName,
			DCM_NameOfPhysiciansReadingStudy,
			DCM_OperatorsName,
			DCM_AdmittingDiagnosesDescription,
			DCM_ReferencedSOPInstanceUID,
			DCM_DerivationDescription,
			DCM_PatientName,
			DCM_PatientID,
			DCM_PatientBirthDate,
			DCM_PatientSex,
			DCM_OtherPatientIDs,
			DCM_OtherPatientNames,
			DCM_PatientAge,
			DCM_PatientSize,
			DCM_PatientWeight,
			DCM_MedicalRecordLocator,
			DCM_EthnicGroup,
			DCM_Occupation,
			DCM_AdditionalPatientHistory,
			DCM_PatientComments,
			DCM_DeviceSerialNumber,
			DCM_ProtocolName,
			DCM_StudyInstanceUID,
			DCM_SeriesInstanceUID,
			DCM_StudyID,
			DCM_FrameOfReferenceUID,
			DCM_SynchronizationFrameOfReferenceUID,
			DCM_ImageComments,
			DCM_RequestAttributesSequence,
			DCM_UID,
			DCM_ContentSequence,
			DCM_StorageMediaFileSetUID,
			DCM_ReferencedFrameOfReferenceUID,
			DCM_RelatedFrameOfReferenceUID,
			DCM_CompressionForce,
		};
		
		for(unsigned int i = 0; i < sizeof(remove_these)/sizeof(*remove_these); ++i){
			remove_all(remove_these[i]);
		}
		
	}
		
	
private:
	void set_pixel_data(uint8_t const* begin, uint8_t const* end);
	void set_pixel_data(uint16_t const* begin, uint16_t const* end);
	
	double get_float64(DcmTagKey key, unsigned int index = 0) const {
		double x;
		assert(sizeof(double) >= 8);
		if(header.getDataset()->findAndGetFloat64(key, x, index).bad()){
			throw mediocr::dicom_key_error("Failed to get some key", key);
		}
		
		return x;
	}
	void set_float64(DcmTagKey key, double value, unsigned int index = 0) {
		assert(sizeof(double) >= 8);
		if(header.getDataset()->putAndInsertFloat64(key, value, index).bad()) {
			throw mediocr::dicom_key_error("Failed to set some key", key);
		}
	}
	float get_float32(DcmTagKey key, unsigned int index = 0) const {
		float x;
		assert(sizeof(float) >= 4);
		if(header.getDataset()->findAndGetFloat32(key, x, index).bad()){
			throw mediocr::dicom_key_error("Failed to get some key", key);
		}
		
		return x;
	}
	void set_float32(DcmTagKey key, float value, unsigned int index = 0) {
		assert(sizeof(float) >= 4);
		if(header.getDataset()->putAndInsertFloat32(key, value, index).bad()) {
			throw mediocr::dicom_key_error("Failed to set some key", key);
		}
	}
	
	std::string get_string(DcmTagKey key) const {
		char const* x;
		if(header.getDataset()->findAndGetString(key, x).bad()){
			throw mediocr::dicom_key_error("Failed to get some key", key);
		}
		
		if(x != NULL){
			return std::string(x);
		} else {
			return std::string("");
		}
	}
	void set_string(DcmTagKey key, std::string value) {
		if(header.getDataset()->putAndInsertString(key, value.c_str()).bad()) {
			throw mediocr::dicom_key_error("Failed to set some key", key);
		}
	}
	
	uint16_t get_uint16(DcmTagKey key) const {
		uint16_t x;
		if(header.getDataset()->findAndGetUint16(key, x).bad()) {
			throw mediocr::dicom_key_error("Failed to get some key", key);
		}
		
		return x;
	}
	void set_uint16(DcmTagKey key, uint16_t value) {
		if(header.getDataset()->putAndInsertUint16(key, value).bad()) {
			throw mediocr::dicom_key_error("Failed to set some key", key);
		}
	}
	
	uint8_t const* get_uint8_array(DcmTagKey key) const {
		uint8_t const* x;
		if(header.getDataset()->findAndGetUint8Array(key, x).bad()){
			throw mediocr::dicom_key_error("Failed to get some key", key);
		}
		
		return x;
	}
	void set_uint8_array(DcmTagKey key, uint8_t const *values, size_t count) {
		if(header.getDataset()->putAndInsertUint8Array(key, values, count).bad()) {
			throw mediocr::dicom_key_error("Failed to set some key", key);
		}
	}
	
	uint16_t const* get_uint16_array(DcmTagKey key) const {
		uint16_t const* x;
		if(header.getDataset()->findAndGetUint16Array(key, x).bad()){
			throw mediocr::dicom_key_error("Failed to get some key", key);
		}
		
		return x;
	}
	void set_uint16_array(DcmTagKey key, uint16_t const *values, size_t count) {
		if(header.getDataset()->putAndInsertUint16Array(key, values, count).bad()) {
			throw mediocr::dicom_key_error("Failed to set some key", key);
		}
	}

	void set_decimal_string(DcmTagKey key, std::vector<float> values){
		std::ostringstream decimal_string;
		auto it = values.begin();
		decimal_string << std::fixed << std::setprecision(3) << *it++;
		for(; it < values.end(); ++it){
			decimal_string << "\\" << std::fixed << std::setprecision(3) << *it;
		}

		set_string(key, decimal_string.str());
	}
	
	mutable DcmFileFormat header;
};

template <>
uint8_t const* dicom_file::get_pixel_data_as<uint8_t const>() const;

template <>
uint16_t const* dicom_file::get_pixel_data_as<uint16_t const>() const;

}
