#include <mediocr-dicom/dicom_file.hpp>
#include <mediocr-dicom/dicom_coders.hpp>

#include <dcmtk/dcmjpeg/djrplol.h>

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
		// Selection Value 7 is chosen because this proved to be the best choice for mammography images
		// (R. Visser, Lossless Compression of Digital Mammograms, IWDM 2006)
		write_to_file(filename, dicom::create_transfer_syntax(EXS_JPEGProcess14TransferSyntax, DJ_RPLossless(7)));
	}

	void dicom_file::write_to_file(std::string filename, E_TransferSyntax syntax, DcmRepresentationParameter* parameters) {
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

		header.getDataset()->chooseRepresentation(syntax, parameters);

		auto const result = header.saveFile(filename.c_str(), syntax);
		if(result.bad()){
			std::string msg{std::string("Error when saving: ") + result.text()};
			throw std::runtime_error(msg);
		}
	}
}
