#include <vector>
#include <iostream>
#include <unistd.h>
#include <complex>
#include <GL/gl.h>
#include <GLFW/glfw3.h>

#define COUNT 1920
#define FREQ 20

#define RENDER_WIDTH 1920
#define RENDER_HEIGHT 1080

void dft(std::vector<std::complex<double>> &in, std::vector<std::complex<double>> &out)
{
	std::fill(out.begin(), out.end(), std::complex<double>(0, 0));
	for (size_t i = 0; i < in.size(); ++i)
	{
		double base = 2 * M_PI * i;
		for (size_t j = 0; j < in.size(); ++j)
		{
			double tmp = base * j / in.size();
			out[i].real(out[i].real() + in[j].real() * cos(tmp));
			out[i].imag(out[i].imag() + in[j].real() * sin(tmp));
		}
	}
}

void dft_re(std::vector<std::complex<double>> &in, std::vector<std::complex<double>> &out)
{
	std::fill(out.begin(), out.end(), std::complex<double>(0, 0));
	for (size_t i = 0; i < in.size(); ++i)
	{
		out[i].real(in[0].real() * cos(0) - in[0].imag() * sin(0));
		double base = 2 * M_PI * i;
		for (size_t j = 1; j < in.size() / 2; ++j)
		{
			double tmp = base * j / in.size();
			out[i].real(out[i].real() + in[j].real() * 2 * cos(tmp) + in[j].imag() * 2 * sin(tmp));
		}
		{
			double tmp = M_PI * i;
			out[i].real(out[i].real() + in[in.size() / 2].real() * cos(tmp) - in[in.size() / 2].imag() * sin(tmp));
		}
		out[i].real(out[i].real() / in.size());
	}
}

void draw(std::vector<std::complex<double>> &in, std::vector<std::complex<double>> &out, std::vector<std::complex<double>> &re)
{
	glBegin(GL_POINTS);
	long prev;
	long value;
	glColor4f(1, 0, 0, 1);
	for (size_t i = 0; i < RENDER_WIDTH; ++i)
	{
		value = RENDER_HEIGHT / 2 + RENDER_HEIGHT / 8 * in[i].real();
		value = std::min((long)RENDER_HEIGHT, std::max((long)0, value));
		glVertex2f(i, value);
		if (i)
		{
			if (value > prev)
			{
				for (int j = prev; j <= value; ++j)
					glVertex2f(i, j);
			}
			else
			{
				for (int j = value; j <= prev; ++j)
					glVertex2f(i, j);
			}
		}
		prev = value;
	}
	glColor4f(0, 1, 0, 1);
	for (size_t i = 0; i < RENDER_WIDTH; ++i)
	{
		value = RENDER_HEIGHT / 2 + RENDER_HEIGHT / 8 * re[i].real();
		value = std::min((long)RENDER_HEIGHT, std::max((long)0, value));
		glVertex2f(i, value);
		if (i)
		{
			if (value > prev)
			{
				for (int j = prev; j <= value; ++j)
					glVertex2f(i, j);
			}
			else
			{
				for (int j = value; j <= prev; ++j)
					glVertex2f(i, j);
			}
		}
		prev = value;
	}
	glColor4f(0, 1, 1, 1);
	for (size_t i = 0; i < RENDER_WIDTH; ++i)
	{
		std::complex<double> &tmp = out[i];
		value = 1 * sqrt(tmp.real() * tmp.real() + tmp.imag() * tmp.imag());
		value = std::min((long)RENDER_HEIGHT, std::max((long)0, value));
		glVertex2f(i, value);
		if (i)
		{
			if (value > prev)
			{
				for (int j = prev; j <= value; ++j)
					glVertex2f(i, j);
			}
			else
			{
				for (int j = value; j <= prev; ++j)
					glVertex2f(i, j);
			}
		}
		prev = value;
	}
	glEnd();
}

int main()
{
	srand(time(nullptr));
	std::vector<std::complex<double>> datas(COUNT);
	std::vector<std::complex<double>> out(COUNT);
	std::vector<std::complex<double>> re(COUNT);
	for (size_t i = 0; i < datas.size(); ++i)
	{
		datas[i] = std::complex<double>(0, 0);
		/*datas[i].real(datas[i].real() + sin(i / static_cast<double>(datas.size()) * FREQ * 1 * M_PI * 2));
		datas[i].real(datas[i].real() + sin(i / static_cast<double>(datas.size()) * FREQ * 2 * M_PI * 2));
		datas[i].real(datas[i].real() + sin(i / static_cast<double>(datas.size()) * FREQ * 3 * M_PI * 2));
		datas[i].real(datas[i].real() + sin(i / static_cast<double>(datas.size()) * FREQ * 4 * M_PI * 2));
		datas[i].real(datas[i].real() + sin(i / static_cast<double>(datas.size()) * FREQ * 5 * M_PI * 2));
		datas[i].real(datas[i].real() + sin(i / static_cast<double>(datas.size()) * FREQ * 6 * M_PI * 2));*/
		//datas[i].real(datas[i].real() + (i % 300) / 300.);
		datas[i].real(datas[i].real() + (i % 1500 < 750 ? 1 : -1));
		//datas[i].real(datas[i].real() + rand() / (float)RAND_MAX);
	}
	dft(datas, out);
	dft_re(out, re);
	glfwInit();
	GLFWwindow *window = glfwCreateWindow(RENDER_WIDTH, RENDER_HEIGHT, "osef", NULL, NULL);
	if (!window)
		return 0;
	glfwMakeContextCurrent(window);
	glEnable(GL_TEXTURE_2D);
	glEnable(GL_BLEND);
	glViewport(0, 0, RENDER_WIDTH, RENDER_HEIGHT);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(0, RENDER_WIDTH, RENDER_HEIGHT, 0, 0, 1);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glClearColor(.1, .1, .1, 1);
	while (!glfwWindowShouldClose(window))
	{
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		glfwWaitEvents();
		draw(datas, out, re);
		glfwSwapBuffers(window);
	}
	return EXIT_SUCCESS;
}
