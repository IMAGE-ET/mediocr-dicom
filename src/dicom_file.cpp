#include <mediocr-dicom/dicom_file.hpp>
#include <mediocr-dicom/dicom_coders.hpp>

#include <dcmtk/dcmjpls/djrparam.h>

namespace mediocr {
	template <>
	uint8_t const* dicom_file::get_pixel_data_as<uint8_t const>() const {
		assert(get_bytes_per_pixel() == 1);
		
		return get_uint8_array(DCM_PixelData);
	}

	template<>
	uint16_t const* dicom_file::get_pixel_data_as<uint16_t const>() const{
		assert(get_bytes_per_pixel() == 2);
		
		return get_uint16_array(DCM_PixelData);
	}

	void dicom_file::set_pixel_data(uint8_t const* begin, uint8_t const* end){
		set_uint8_array(DCM_PixelData, begin, std::distance(begin, end));
	}

	void dicom_file::set_pixel_data(uint16_t const* begin, uint16_t const* end){
		set_uint16_array(DCM_PixelData, begin, std::distance(begin, end));
	}

	void dicom_file::write_to_file(std::string filename) {
		// JPEG-LS Lossless gets really nice compression and it's lossless!
		write_to_file(filename, dicom::create_transfer_syntax(EXS_JPEGLSLossless, DJLSRepresentationParameter(0, true)));
	}

	void dicom_file::write_to_file(std::string filename, E_TransferSyntax syntax, DcmRepresentationParameter* parameters) {
		update_content_date_time();

		auto result = header.getDataset()->chooseRepresentation(syntax, parameters);
		if(result.bad()) {
			std::string msg{std::string("Error when converting to representation: ") + result.text()};
			throw std::runtime_error(msg);
		}

		result = header.saveFile(filename.c_str(), syntax);
		if(result.bad()){
			std::string msg{std::string("Error when saving: ") + result.text()};
			throw std::runtime_error(msg);
		}
	}

	void dicom_file::generate_series_instance_uid(const char* site_root){
		char uid[1024];
		dcmGenerateUniqueIdentifier(uid, site_root);
		set_series_uid(uid);
	}

	void dicom_file::generate_study_instance_uid(const char* site_root){
		char uid[1024];
		dcmGenerateUniqueIdentifier(uid, site_root);
		set_study_uid(uid);
	}

	void dicom_file::generate_sop_instance_uid(const char* site_root){
		char uid[1024];
		dcmGenerateUniqueIdentifier(uid, site_root);
		set_instance_uid(uid);
	}

	void dicom_file::update_content_date_time(){
		time_t rawtime;
		time(&rawtime);
		struct tm* timeinfo = localtime(&rawtime);

		char tm[7] = {0};
		sprintf(tm, "%02d%02d%02d", timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec);

		// the tm structure defines month as 0-11
		// the tm structure defines mday as 1-31
		char dt[9] = {0};
		sprintf(dt, "%04d%02d%02d", timeinfo->tm_year + 1900, timeinfo->tm_mon + 1, timeinfo->tm_mday);

		set_string(DCM_ContentDate, dt);
		set_string(DCM_ContentTime, tm);
	}
}
