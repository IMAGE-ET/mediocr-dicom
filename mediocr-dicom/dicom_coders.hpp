#pragma once

#include <string>

namespace mediocr {

/*! Returns true if the file has the Dicom magic */
bool file_has_dicom_magic(std::string filepath);

namespace dicom {
	struct Coders {
		Coders();
		~Coders();

	private:
		static void djencoder_register();
		static void djencoder_cleanup();

		static void djdecoder_register();
		static void djdecoder_cleanup();
	};
}

}
