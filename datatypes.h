#include <iio.h>

struct extra_info {
	struct iio_device *dev;
	float *data_ref;
	off_t offset;
	int shadow_of_enabled;
	bool may_be_enabled;
	double lo_freq;
};

struct extra_dev_info {
	bool input_device;
	struct iio_buffer *buffer;
	unsigned int sample_count;
	unsigned int buffer_size;
	unsigned int channel_trigger;
	bool channel_trigger_enabled;
	bool trigger_falling_edge;
	float trigger_value;
	double adc_freq;
	char adc_scale;
	float **channels_data_copy;
	//GSList *plots_sample_counts;
	float plugin_fft_corr;
};
