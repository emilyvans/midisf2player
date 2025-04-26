#ifndef INCLUDE_MIDI_FILE_HPP_
#define INCLUDE_MIDI_FILE_HPP_
#include <optional>
#include <stdint.h>
#include <variant>
#include <vector>

struct MIDI_DIVISION {
	bool format;
	union {
		uint16_t ticks_per_quarter_note;
		struct {
			uint8_t frames_per_second;
			uint8_t ticks_per_frame;
		};
	};
};

struct MIDI_HEADER {
	char type[5] = {0};
	uint32_t length;
	uint16_t format;
	uint16_t ntracks;
	MIDI_DIVISION division;
};

struct MIDI_SYSEX_EVENT {
	std::vector<uint8_t> data;
};

struct MIDI_META_EVENT {
	uint8_t type;
	std::vector<uint8_t> data;
};

struct MIDI_MIDI_EVENT {
	uint8_t type;
	uint8_t channel;
	std::vector<uint8_t> data;
};

struct MIDI_EVENT {
	uint32_t delta_time;
	uint8_t type;
	std::variant<MIDI_SYSEX_EVENT, MIDI_META_EVENT, MIDI_MIDI_EVENT> event;
};

struct MIDI_TRACK {
	char type[5] = {0};
	uint32_t length;
	std::vector<MIDI_EVENT> events;
};

struct MIDI_FILE {
	MIDI_HEADER header;
	std::vector<MIDI_TRACK> tracks;
};

std::optional<MIDI_FILE> read_midi_file(const char *filename);

#endif // INCLUDE_MIDI_FILE_HPP_
