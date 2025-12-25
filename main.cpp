#include <iostream>
#include <iomanip>
#include <vector>
#include <complex>
#include <cmath>
#include <chrono> 

const double PI = 3.14159265358979323846;

const int POWER_N = 11;             
const int N = 1 << POWER_N;         
const double A = 2.49;
const double B = 0.16;
const double W1 = 3.0;
const double W2 = 199.0;
const double PHI = PI / 4.0;        

using Complex = std::complex<double>;

void DFT(const std::vector<Complex>& Data, std::vector<Complex>& Result)
{
    int size = Data.size();
    Result.assign(size, { 0, 0 });
    for (int m = 0; m < size; m++) {
        for (int n = 0; n < size; n++) {
            double angle = -2.0 * PI * m * n / size;
            Result[m] += Data[n] * std::polar(1.0, angle);
        }
    }
}

void IDFT(const std::vector<Complex>& Data, std::vector<Complex>& Result)
{
    int size = Data.size();
    Result.assign(size, { 0, 0 });
    for (int n = 0; n < size; n++) {
        for (int m = 0; m < size; m++) {
            double angle = 2.0 * PI * m * n / size;
            Result[n] += Data[m] * std::polar(1.0, angle);
        }
        Result[n] /= (double)size;
    }
}

void FFT(const std::vector<Complex>& Data, std::vector<Complex>& Result)
{
    int size = Data.size();
    int M = size / 2;
    Result.clear(); Result.resize(size);
    Complex Exp, U, V;

    for (int m = 0; m < M; m++) {
        U = { 0.0, 0.0 };
        V = { 0.0, 0.0 };
        for (int n = 0; n < M; n++) {
            double angle1 = -2.0 * PI * m * n / M;
            Exp = std::polar(1.0, angle1);
            U += Data[2 * n] * Exp;
            V += Data[2 * n + 1] * Exp;
        }
        double angle2 = -2.0 * PI * m / size;
        Exp = std::polar(1.0, angle2);
        Result[m] = U + Exp * V;
        Result[m + M] = U - Exp * V;
    }
}

void IFFT(const std::vector<Complex>& Data, std::vector<Complex>& Result)
{
    int size = Data.size();
    std::vector<Complex> Conj(size);
    for (int i = 0; i < size; i++) Conj[i] = std::conj(Data[i]);
    std::vector<Complex> Temp;
    FFT(Conj, Temp);
    Result.assign(size, 0);
    for (int i = 0; i < size; i++) {
        Result[i] = std::conj(Temp[i]) / (double)size;
    }
}

int main()
{
    setlocale(LC_ALL, "Russian");
    std::vector<Complex> Z(N);
    std::vector<Complex> DFT_Result(N);
    std::vector<Complex> FFT_Result(N);
    std::vector<Complex> Filtered_Z(N);

    for (int j = 0; j < N; j++) {
        double arg1 = 2.0 * PI * W1 * j / N + PHI;
        double arg2 = 2.0 * PI * W2 * j / N;
        double val = A * cos(arg1) + B * cos(arg2);
        Z[j] = Complex(val, 0.0);
    }

    auto start = std::chrono::high_resolution_clock::now();
    DFT(Z, DFT_Result);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff_dft = end - start;
    std::cout << "Время DFT: " << diff_dft.count() << " s" << std::endl;

    start = std::chrono::high_resolution_clock::now();
    FFT(Z, FFT_Result);
    end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff_fft = end - start;
    std::cout << "Время FFT: " << diff_fft.count() << " s" << std::endl;
    std::cout << "----------------------------------------------------" << std::endl;

    int w = 15;
    std::cout << std::left << std::setw(6) << "m" << std::setw(w) << "Re z" << std::setw(w) << "Re z^"
        << std::setw(w) << "Im z^" << std::setw(w) << "Амплитуда |z^|" << std::setw(w) << "Фаза phi" << std::endl;

    double threshold = 0.001;
    for (int m = 0; m < N; m++) {
        double amplitude = std::abs(FFT_Result[m]);
        if (amplitude > threshold) {
            std::cout << std::left << std::setw(6) << m << std::setw(w) << Z[m].real() << std::setw(w) << FFT_Result[m].real()
                << std::setw(w) << FFT_Result[m].imag() << std::setw(w) << amplitude << std::setw(w) << std::arg(FFT_Result[m]) << std::endl;
        }
    }
    std::cout << "----------------------------------------------------" << std::endl;

    for (int m = 0; m < N; m++) {
        if (std::abs(m - W2) < 0.5 || std::abs(m - (N - W2)) < 0.5) {
            FFT_Result[m] = { 0.0, 0.0 };
        }
    }

    std::cout << "\nПосле фильтрации" << std::endl;
    std::cout << std::left << std::setw(6) << "m" << std::setw(w) << "Rez^" << std::setw(w) << "Амплитуда |z^|" << std::endl;
    for (int m = 0; m < N; m++) {
        double amplitude = std::abs(FFT_Result[m]);
        if (amplitude > threshold) {
            std::cout << std::left << std::setw(6) << m << std::setw(w) << FFT_Result[m].real() << std::setw(w) << amplitude << std::endl;
        }
    }
    std::cout << "----------------------------------------------------" << std::endl;

    IFFT(FFT_Result, Filtered_Z);
    std::cout << "Сравнение" << std::endl;
    std::cout << std::left << std::setw(6) << "j" << std::setw(w) << "Ориг Rez" << std::setw(w) << "Фильтр Rez" << std::setw(w) << "Разница" << std::endl;
    for (int j = 0; j < 10; j++) {
        std::cout << std::left << std::setw(6) << j << std::setw(w) << Z[j].real() << std::setw(w) << Filtered_Z[j].real() << std::setw(w) << (Z[j].real() - Filtered_Z[j].real()) << std::endl;
    }
    std::cout << "----------------------------------------------------" << std::endl;

    std::vector<Complex> Z_piecewise(N), FFT_piece(N);
    for (int j = 0; j < N; j++) {
        if ((j >= N / 4 && j <= N / 2) || (j > 3 * N / 4 && j < N)) {
            Z_piecewise[j] = Complex(A + B * cos(2.0 * PI * W2 * j / N), 0.0);
        }
        else {
            Z_piecewise[j] = Complex(0.0, 0.0);
        }
    }
    FFT(Z_piecewise, FFT_piece);

    std::cout << "Спектр" << std::endl;
    std::cout << std::left << std::setw(6) << "m" << "Амплитуда |z^|" << std::endl;
    for (int m = 0; m < 30; m++) {
        std::cout << std::left << std::setw(6) << m << std::abs(FFT_piece[m]) << std::endl;
    }
    std::cout << "----------------------------------------------------" << std::endl;

    return 0;
}