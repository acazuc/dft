#include <vector>
#include <iostream>
#include <unistd.h>
#include <complex>
#include <GL/gl.h>
#include <GLFW/glfw3.h>

#define COUNT 1920
#define FREQ 2

#define RENDER_WIDTH 1920
#define RENDER_HEIGHT 1080

std::vector<std::complex<double>> g_in(COUNT);
std::vector<std::complex<double>> g_out(COUNT);
std::vector<std::complex<double>> g_re(COUNT);
float g_scale = 8;

template <typename T>
void dft(std::vector<std::complex<T>> &in, std::vector<std::complex<T>> &out)
{
	out.resize(in.size());
	std::fill(out.begin(), out.end(), std::complex<T>(0, 0));
	for (size_t i = 0; i < out.size(); ++i)
	{
		T base = 2 * M_PI * i / in.size();
		T tmp = 0;
		for (size_t j = 0; j < in.size(); ++j)
		{
			out[i].real(out[i].real() + in[j].real() * cos(tmp));
			out[i].imag(out[i].imag() + in[j].real() * sin(tmp));
			tmp += base;
		}
	}
}

template <typename T>
void idft(std::vector<std::complex<T>> &in, std::vector<std::complex<T>> &out)
{
	out.resize(in.size());
	std::fill(out.begin(), out.end(), std::complex<T>(0, 0));
	for (size_t i = 0; i < out.size(); ++i)
	{
		out[i].real(in[0].real() * cos(0) + in[0].imag() * sin(0));
		T base = 2 * M_PI * i / in.size();
		T tmp = 0;
		for (size_t j = 1; j < in.size() / 2; ++j)
		{
			tmp += base;
			out[i].real(out[i].real() + 2 * (in[j].real() * cos(tmp) + in[j].imag() * sin(tmp)));
		}
		{
			tmp = M_PI * i;
			out[i].real(out[i].real() + in[in.size() / 2].real() * cos(tmp) + in[in.size() / 2].imag() * sin(tmp));
		}
		out[i].real(out[i].real() / in.size());
	}
}

template <typename T>
void dct(std::vector<std::complex<T>> &in, std::vector<std::complex<T>> &out)
{
	out.resize(in.size());
	std::fill(out.begin(), out.end(), std::complex<T>(0, 0));
	T c1 = M_PI / out.size();
	for (size_t i = 0; i < out.size(); ++i)
	{
		for (size_t j = 0; j < in.size(); ++j)
		{
			out[i].real(out[i].real() + in[j].real() * cos(c1 * (j + .5) * i));
		}
	}
}

template <typename T>
void dct(std::vector<T> &in, std::vector<T> &out)
{
	out.resize(in.size());
	std::fill(out.begin(), out.end(), 0);
	T c1 = M_PI / out.size();
	for (size_t i = 0; i < out.size(); ++i)
	{
		for (size_t j = 0; j < in.size(); ++j)
		{
			out[i] += in[j] * cos(c1 * (j + .5) * i);
		}
	}
}

template <typename T>
void idct(std::vector<std::complex<T>> &in, std::vector<std::complex<T>> &out)
{
	out.resize(in.size());
	T c1 = M_PI / in.size();
	for (size_t i = 0; i < out.size(); ++i)
	{
		out[i].real(.5 * in[0].real());
		for (size_t j = 1; j < in.size(); ++j)
		{
			out[i].real(out[i].real() + in[j].real() * cos(c1 * j * (i + .5)));
		}
	}
}

template <typename T>
void idct(std::vector<T> &in, std::vector<T> &out)
{
	out.resize(in.size());
	T c1 = M_PI / in.size();
	for (size_t i = 0; i < out.size(); ++i)
	{
		out[i] = .5 * in[0];
		for (size_t j = 1; j < in.size(); ++j)
		{
			out[i] += in[j] * cos(c1 * j * (i + .5));
		}
	}
}

template <typename T>
void mdct(std::vector<T> &in, std::vector<T> &out)
{
	out.resize(in.size() / 2);
	std::fill(out.begin(), out.end(), 0);
	T c1 = M_PI / out.size();
	T c2 = .5 + out.size() / 2;
	for (size_t i = 0; i < out.size(); ++i)
	{
		for (size_t j = 0; j < in.size(); ++j)
		{
			out[i] += in[j] * cos(c1 * (j + c2) * (i + .5));
		}
	}
}

template <typename T>
void mdct(std::vector<std::complex<T>> &in, std::vector<std::complex<T>> &out)
{
	out.resize(in.size() / 2);
	std::fill(out.begin(), out.end(), std::complex<T>(0, 0));
	T c1 = M_PI / out.size();
	T c2 = .5 + out.size() / 2;
	for (size_t i = 0; i < out.size(); ++i)
	{
		for (size_t j = 0; j < in.size(); ++j)
		{
			out[i].real(out[i].real() + in[j].real() * cos(c1 * (j + c2) * (i + .5)));
		}
	}
}

template <typename T>
void imdct(std::vector<T> &in, std::vector<T> &out)
{
	out.resize(in.size() * 2);
	std::fill(out.begin(), out.end(), 0);
	T c1 = M_PI / in.size();
	T c2 = .5 + in.size() / 2;
	for (size_t i = 0; i < out.size(); ++i)
	{
		for (size_t j = 0; j < in.size(); ++j)
		{
			out[i] += in[j] * cos(c1 * (i + c2) * (j + .5));
		}
		out[i] /= in.size();
	}
}

template <typename T>
void imdct(std::vector<std::complex<T>> &in, std::vector<std::complex<T>> &out)
{
	out.resize(in.size() * 2);
	std::fill(out.begin(), out.end(), std::complex<T>(0, 0));
	T c1 = M_PI / in.size();
	T c2 = .5 + in.size() / 2;
	for (size_t i = 0; i < out.size(); ++i)
	{
		for (size_t j = 0; j < in.size(); ++j)
		{
			out[i].real(out[i].real() + in[j].real() * cos(c1 * (i + c2) * (j + .5)));
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
		value = RENDER_HEIGHT / 2 + RENDER_HEIGHT / g_scale * in[i].real();
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
		value = RENDER_HEIGHT / 2 + RENDER_HEIGHT / g_scale * re[i].real();// / 960;
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

void run()
{
	dft(g_in, g_out);
	idft(g_out, g_re);
}

static bool g_changed = false;
int last_y = -1;
int last_x = -1;

void cursor_pos_callback(GLFWwindow* window, double xpos, double ypos)
{
	int state = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT);
	if (state != GLFW_PRESS)
		return;
	int x = xpos;
	if (x < 0 || x > RENDER_WIDTH)
		return;
	if (x == last_x)
		return;
	int y = ypos;
	if (y < 0 || y > RENDER_HEIGHT)
		return;
	if (last_x == -1)
		g_in[x] = (y + 20 - (RENDER_HEIGHT / 2)) / RENDER_HEIGHT * g_scale;
	else
	{
		if (x > last_x)
		{
			for (int i = last_x; i <= x; ++i)
				g_in[i] = (last_y + (y - last_y) * ((i - last_x) / (float)(x - last_x)) + 20 - (RENDER_HEIGHT / 2)) / RENDER_HEIGHT * g_scale;
		}
		else
		{
			for (int i = x; i <= last_x; ++i)
				g_in[i] = (y + (last_y - y) * ((i - x) / (float)(last_x - x)) + 20 - (RENDER_HEIGHT / 2)) / RENDER_HEIGHT * g_scale;
		}
	}
	last_x = x;
	last_y = y;
	g_changed = true;
}

void mouse_button_callback(GLFWwindow* window, int button, int action, int mods)
{
	if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_RELEASE && g_changed)
		run();
	last_x = -1;
}

int main()
{
	srand(time(nullptr));
	for (size_t i = 0; i < g_in.size(); ++i)
	{
		g_in[i] = std::complex<double>(0, 0);
		int osef = 2;
		for (int j = 0; j < osef; ++j)
			g_in[i].real(g_in[i].real() + sin(i / static_cast<double>(g_in.size()) * FREQ * pow(2, j) * M_PI * 2) / osef);
		//g_in[i].real(g_in[i].real() + (i % 300) / 300.);
		//g_in[i].real(g_in[i].real() + (i < RENDER_WIDTH / 2 ? 1 : -1));
		//g_in[i].real(g_in[i].real() + rand() / (float)RAND_MAX);
	}
	run();
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
	glfwSetCursorPosCallback(window, cursor_pos_callback);
	glfwSetMouseButtonCallback(window, mouse_button_callback);
	while (!glfwWindowShouldClose(window))
	{
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		glfwWaitEvents();
		draw(g_in, g_out, g_re);
		glfwSwapBuffers(window);
	}
	return EXIT_SUCCESS;
}
