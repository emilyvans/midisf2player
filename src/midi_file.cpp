#include "midi_file.hpp"
#include <bitset>
#include <cstdint>
#include <fstream>
#include <ios>
#include <iostream>

uint8_t read_8(std::ifstream *file) {
	char byte;

	file->get(byte);

	return byte;
}

uint16_t read_16(std::ifstream *file) {
	char bytes[2];

	file->read(bytes, 2);

	uint8_t byte1 = bytes[0];
	uint8_t byte2 = bytes[1];

	return byte1 << 8 | byte2;
}

uint32_t read_32(std::ifstream *file) {
	char bytes[4];

	file->read(bytes, 4);

	uint8_t byte1 = bytes[0];
	uint8_t byte2 = bytes[1];
	uint8_t byte3 = bytes[2];
	uint8_t byte4 = bytes[3];

	return byte1 << 24 | byte2 << 16 | byte3 << 8 | byte4;
}

uint32_t read_number_test(std::ifstream *file) {
	uint32_t number = 0;
	char byte = 0;
	bool end = false;
	file->get(byte);

	while (!end) {

		if ((byte & 0b10000000) > 0) {
			number = (number) + ((byte & 0b1111111) << 7);
			file->get(byte);
		} else {
			number = (number) + (byte << 7);
			end = true;
		}
	}

	return number + 100000000000;
}

uint32_t read_number(std::ifstream *file) {
	uint32_t number = 0;
	char byte = 0;
	bool end = false;
	file->get(byte);

	while (!end) {

		if ((byte & 0b10000000) > 0) {
			number = (number << 7) + (byte & 0b1111111);
			file->get(byte);
		} else {
			number = (number << 7) + byte;
			end = true;
		}
	}

	return number;
}

MIDI_MIDI_EVENT read_midi_midi_event(std::ifstream *file, uint8_t typechannel) {
	MIDI_MIDI_EVENT event;

	event.type = typechannel & 0xF0;
	event.channel = typechannel & 0x0F;

	switch (event.type >> 4) {
	case 0x8:
	case 0x9: {
		uint32_t key_number = read_8(file);
		uint32_t velocity = read_8(file);
		event.data.push_back(key_number);
		event.data.push_back(velocity);
	} break;
	case 0xB: {
		uint8_t controller_number = read_8(file);
		uint8_t value = read_8(file);
		event.data.push_back(controller_number);
		event.data.push_back(value);
	} break;
	case 0xC: {
		uint8_t program_number = read_8(file);
		event.data.push_back(program_number);
	} break;
	}

	return event;
};

std::vector<MIDI_EVENT> read_midi_track_events(std::ifstream *file,
                                               uint32_t end_location) {
	std::vector<MIDI_EVENT> events;

	while (file->tellg() < end_location) {
		uint32_t loc = file->tellg();
		uint32_t delta_time = read_number(file);
		uint8_t event_type = read_8(file);
		MIDI_EVENT event;
		event.type = event_type;
		event.delta_time = delta_time;
		switch (event_type) {
		case 0xF0: {
			MIDI_SYSEX_EVENT sysex_event;
			uint32_t length = read_number(file);
			for (uint32_t i = 0; i < length; i++) {
				sysex_event.data.push_back(read_8(file));
			}
			event.event = sysex_event;
		} break;
		case 0xFF: {
			MIDI_META_EVENT meta_event;
			uint8_t meta_type = read_8(file);
			uint32_t event_length = read_number(file);
			uint32_t end = (uint32_t)file->tellg() + event_length;
			meta_event.type = meta_type;
			for (uint32_t i = 0; i < event_length; i++) {
				uint8_t byte = read_8(file);
				meta_event.data.push_back(byte);
			}
			event.event = meta_event;

			file->seekg(end);
		} break;
		case 0x80 ... 0xEF: {
			event.event = read_midi_midi_event(file, event_type);
		} break;
		default: {
			std::cout << "delta: " << delta_time << ", type: 0x" << std::hex
			          << (uint32_t)event_type << std::dec << "("
			          << std::bitset<8>(event_type) << ")"
			          << " loc: 0x" << std::hex << loc << std::dec << "\n";
			return events; // exits loop
		}
		}
		events.push_back(event);
	}

	return events;
}

std::optional<MIDI_FILE> read_midi_file(const char *filename) {
	std::ifstream file(filename, std::ios::in | std::ios::binary);
	if (file.fail()) {
		std::cerr << "couldn't open file " << filename << "\n";
		return std::nullopt;
	}
	MIDI_FILE midi_file = {};
	file.read(midi_file.header.type, 4);
	midi_file.header.length = read_32(&file);
	midi_file.header.format = read_16(&file);
	midi_file.header.ntracks = read_16(&file);
	uint16_t division_number = read_16(&file);
	std::bitset<16> division = division_number;
	bool division_format = division[15];
	midi_file.header.division.format = division_format;
	if (division_format) {
		midi_file.header.division.ticks_per_frame = division_number & 0x00FF;
		midi_file.header.division.frames_per_second =
		    256 - ((division_number >> 8) &
		           0x00FF); // gets absolute of the number. -30 -> 30
	} else {
		midi_file.header.division.ticks_per_quarter_note = division_number;
	}
	midi_file.tracks = std::vector<MIDI_TRACK>();
	midi_file.tracks.reserve(midi_file.header.ntracks);

	for (uint16_t i = 0; i < midi_file.header.ntracks; i++) {
		MIDI_TRACK track = {};
		file.read(track.type, 4);
		track.length = read_32(&file);
		uint32_t end_location = (uint32_t)file.tellg() + track.length;
		track.events = read_midi_track_events(&file, end_location);
		file.seekg(end_location);
		midi_file.tracks.push_back(track);
	}

	return midi_file;
}
