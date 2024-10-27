#ifndef SF2_FILE_HPP_
#define SF2_FILE_HPP_
#include <cstdint>
#include <cstdlib>
#include <cwchar>
#include <string>
#include <vector>

struct sf_preset_header { // sfPresetHeader
	char preset_name[20];
	uint16_t preset;
	uint16_t bank;
	uint16_t perset_bag_ndx;
	uint32_t library;
	uint32_t genre;
	uint32_t morphology;
};

struct sf_preset_zone { // sfPresetBag
	uint16_t GenNdx;
	uint16_t ModNdx;
};

struct sf_mod_list { // sfModList
	SF_modulator src_oper;
	SF_generator dest_oper;
	uint16_t mod_amount;
	SF_modulator amt_src_oper;
	SF_transform trans_oper;
};

class SF2_FILE {
  public:
	SF2_FILE(std::string file_path);

  private:
  public:
  private:
	std::string file_path;
	uint32_t file_length;
	// INFO
	uint32_t major_version;
	uint32_t minor_version;
	std::string sound_engine;
	std::string name;
	std::string date_of_creation;
	std::string tool;
	// sdta
	std::vector<uint8_t> samples;
	// pdta
	std::vector<sf_preset_header> preset_headers;
};

#endif // SF2_FILE_HPP
