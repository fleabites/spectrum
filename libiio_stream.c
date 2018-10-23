
#include <iio.h>
#include <math.h>
#include <limits.h>

float dither(float f)
{
	float r1 = (float)rand() / (float)RAND_MAX;
	float r2 = (float)rand() / (float)RAND_MAX;

	return f + (r1 - r2) * 16.0f;
}

int main(void)
{
	struct iio_context *ctx;
	struct iio_device *dds;
	struct iio_buffer *txbuf;
	unsigned int i, num_channels;
	int16_t *buf, *sine;

	ctx = iio_create_local_context();
	if (!ctx) {
		perror("Could not create IIO context");
		return -1;
	}

	dds = iio_context_find_device(ctx, "cf-ad9361-dds-core-lpc");
	if (!dds) {
		perror("Could not find device");
		return -1;
	}

	num_channels = iio_device_get_channels_count(dds);
	for (i = 0; i < num_channels; i++)
		iio_channel_enable(iio_device_get_channel(dds, i));

	txbuf = iio_device_create_buffer(dds, 1024 * 256, 0);
	if (!txbuf) {
		perror("Could not allocate buffer");
		return -1;
	}

	sine = malloc(sizeof(int16_t) * 1024 * 256 * 2);

	for (i = 0; i < 1024 * 256; i++) {
		sine[i*2 + 0] = dither(cos(2 * M_PI * i / 256.0) * 0x4000);
		sine[i*2 + 1] = dither(sin(2 * M_PI * i / 256.0) * 0x4000);
	}

	while (1) {
		buf = iio_buffer_start(txbuf);
		memcpy(buf, sine, 1024 * 256 * sizeof(int16_t) * 2);
		iio_buffer_push(txbuf);
	}

	return 0;
}
