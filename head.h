#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <math.h>
#include <cstdio>
#include <cstring>
#include <sstream>
#include <iomanip>
#include <stdexcept> // Для runtime_error

using namespace std;

/* Функция для решения системы линейных уравнений методом сопряженных градиентов */
void conjugate_gradient_s_limit (int size, double **A, double *B, double *x_conju, int* index_v, double* react,
                                 int iter_max = 1e4, double epsilon = 1e-14, bool showing_param = 0)
{
  double *p = new double[size];        // p_k направление спуска 
  double *r = new double[size];        // r_k невязка
  double *tt = new double[size];       // tt = A * p_k 
  double alpha, beta, tmp1, nevyzka, B_norm = 0;
  
  // возвращает 1 если позиция i в блок - векторе 2*L в тех индексах, которые надо занулить
  auto is_in_index_v = [](int i, int* index_v) 
     { return (i == 2 * index_v[0] || i == 2 * index_v[1] || i == 2 * index_v[2] || i == 2 * index_v[3] || i == (2 * index_v[3] + 1));}; 
  
  // зануляет позиции опор
  auto cond_vect = [](double* a, int size, int* index_v) { for (size_t i = 0; i < 4; i++)  { a[2 * index_v[i]] = 0.0; } a[2 * index_v[3] + 1] = 0.0; }; 

  cond_vect (x_conju, size, index_v);

  /* Инициализация r_0 и p_0 */
  for (int i = 0; i < size; i++)
    {
      r[i] = B[i];
      for (size_t k = 0; k < size; k++)  { r[i] -= A[i][k] * x_conju[k]; }  // случай ненулевого начального условия x_0
      p[i] = r[i];
      B_norm += B[i] * B[i];  /* Норма B */
    }

  cond_vect (p, size, index_v);
  B_norm = sqrt(B_norm);

  /* Итерации метода сопряженных градиентов */
  for (int count = 0; count < iter_max; count++)  /* Счетчик итераций */
    {
      alpha = 0.0;
      tmp1 = 0.0;
      beta = 0.0;

      /* Вычисление A * p_k */
      for (int i = 0; i < size; i++)
        {
          tt[i] = 0.0;
          if (!is_in_index_v (i, index_v))
            {
              for (int k = 0; k < size; k++)  /* Вычисляем tt = A * p_k */
                {
                  tt[i] += A[i][k] * p[k];
                }
            }
          tmp1 += r[i] * r[i];    /* tmp1 = (r_k, r_k) */
          alpha += p[i] * tt[i];  /* alpha = (p_k, A * p_k) */
        }

      /* Проверка деления на ноль */
      if (fabs(alpha) < 1e-12)  { cout << "Ошибка: деление на ноль в alpha" << endl;    break; }

      alpha = tmp1 / alpha;

      /* Обновление x_k+1 и r_k+1 */
      for (int i = 0; i < size; i++)
        {
          x_conju[i] += alpha * p[i];
          r[i] -= alpha * tt[i];
        }

      /* Вычисление beta */
      for (int i = 0; i < size; i++)
        {
          beta += r[i] * r[i];
        }

      /* Норма остатка */
      nevyzka = beta;
      beta = beta / tmp1;

      /* Метод Флетчера - Ривса | нужен рестарт на каждом size+1 иттерации нужно забывать направление спуска */
      if (count % (size + 1) == 0)
        {
          beta = 0.0;
        }

      /* Обновление p_k+1 */
      for (int i = 0; i < size; i++)
        {
          p[i] = r[i] + beta * p[i];
        }
      cond_vect (p, size, index_v);

      /* Проверка условия останова */
      if (nevyzka < epsilon * B_norm)
        {
          if (showing_param)
            {
              cout << "Метод сопряженных градиентов завершился после " << count + 1 << " итераций." << endl;
            }
          break;
        }
    }
  // after itterations
  for (int i = 0; i < size; i++)
    {
      r[i] = -B[i];
      for (size_t k = 0; k < size; k++)  { r[i] += A[i][k] * x_conju[k]; } 
    }
  //cout << "r = ";   show_vect (r, size, 0);
  
  for (size_t i = 0; i < 4; i++)
    {
      react[i] = r[2 * index_v[i]];
    }
  react[4] = r[2 * index_v[3] + 1];

  delete[] p;
  delete[] r;
  delete[] tt;
}

/* Функция для решения системы линейных уравнений методом сопряженных градиентов */
void conjugate_gradient (int size, double **A, double *B, double *x_conju, int iter_max = 1e4, double epsilon = 1e-14, bool showing_param = 0)
{
  double *p = new double[size];    /* p_k */
  double *r = new double[size];    /* r_k */
  double *tt = new double[size];   /* tt = A * p_k */
  double alpha, beta, tmp1, nevyzka, B_norm = 0;

  /* Инициализация r_0 и p_0 */
  for (int i = 0; i < size; i++)
    {
      r[i] = B[i];
      for (size_t k = 0; k < size; k++)  { r[i] -= A[i][k] * x_conju[k]; }  // случай ненулевого начального условия x_0
      p[i] = r[i];
      B_norm += B[i] * B[i];  /* Норма B */
    }

  B_norm = sqrt(B_norm);

  /* Итерации метода сопряженных градиентов */
  int count = 0;
  for (count = 0; count < iter_max; count++)  /* Счетчик итераций */
    {
      alpha = 0.0;
      tmp1 = 0.0;
      beta = 0.0;

      /* Вычисление A * p_k */
      for (int i = 0; i < size; i++)
        {
          tt[i] = 0.0;
          for (int k = 0; k < size; k++)  /* Вычисляем tt = A * p_k */
            {
              tt[i] += A[i][k] * p[k];
            }

          tmp1 += r[i] * r[i];  /* tmp1 = (r_k, r_k) */
          alpha += p[i] * tt[i];  /* alpha = (p_k, A * p_k) */
        }

            /* Проверка деления на ноль */
      if (fabs(alpha) < 1e-12)
        {
          cout << "Ошибка: деление на ноль в alpha" << endl;
          break;
        }

      /* Вычисление alpha */
      alpha = tmp1 / alpha;

      /* Обновление x_k+1 и r_k+1 */
      for (int i = 0; i < size; i++)
        {
          x_conju[i] += alpha * p[i];
          r[i] -= alpha * tt[i];
        }

      /* Вычисление beta */
      for (int i = 0; i < size; i++)
        {
          beta += r[i] * r[i];
        }

      // cout << "Итерация " << count << " Невязка: " << sqrt(beta) << endl;

      /* Норма остатка */
      nevyzka = beta;
      beta = beta / tmp1;

      /* Еще он звисит но начальной точки */
      /* Метод Флетчера - Ривса нужен рестарт на каждом size+1 иттерации нужно забывать направление спуска */
      if (count % (size + 1) == 0)
        {
          beta = 0.0;
        }

      /* Обновление p_k+1 */
      for (int i = 0; i < size; i++)
        {
          p[i] = r[i] + beta * p[i];
        }

      /* Проверка условия останова */
      if (nevyzka < epsilon * B_norm)
        {
          if (showing_param)
            {
              cout << "Метод сопряженных градиентов завершился после " << count + 1 << " итераций." << endl;
            }
          break;
        }
    }

  delete[] p;
  delete[] r;
  delete[] tt;
}

void LU_(int size, double **A, double *B, double* x_LU)
{
  double** L = new double *[size];                                                                               
  double** U = new double *[size];        

  for (size_t i = 0; i < size; i++)
    {
      U[i] = new double[size];   L[i] = new double[size]; 
    }        

  for (size_t i = 0; i < size; i++) 
   {
      /* Вычисление элементов U для строки i */
      for (size_t j = i; j < size; j++) 
        {
          U[i][j] = A[i][j];
          for (size_t k = 0; k < i; k++) 
            {
              U[i][j] -= L[i][k] * U[k][j];
            }
        }

      /* Вычисление элементов L для столбца i */
      for (size_t j = i + 1; j < size; j++) 
        {
          L[j][i] = A[j][i];
          for (size_t k = 0; k < i; k++) 
            {
              L[j][i] -= L[j][k] * U[k][i];
            }
          L[j][i] /= U[i][i];
        }

      /* Диагональный элемент L[i][i] всегда равен 1 */
      L[i][i] = 1.0;
  }
  
  for (int i = 0; i < size; i++)                   // обратный ход L
    {
      x_LU[i] = B[i];

      for (size_t j = 0; j < i; j++)
        {
          x_LU[i] -= L[i][j] * x_LU[j];
        }
    }

  for (int i = size - 1; i > -1; i--)             // обратный ход U
    {
      for (size_t j = i + 1; j < size; j++)
        {
          x_LU[i] -= U[i][j] * x_LU[j];
        }
      x_LU[i] /= U[i][i];
    }

  for (size_t i = 0; i < size; i++) 
    {
      delete[] L[i];
      delete[] U[i];
    }
  delete[] L;   delete[] U;  
}

/* Вычисление значений базисной функции и её производных */
double base_function (int p, double t, double h, int derivative)
{
  if (t < 0 || t > 1)
    {
      return 0.0;
    }

  switch (derivative)
    {
      case 0: // Базисные функции
        switch (p)
          {
            case 1: return 1 - 3 * pow(t, 2) + 2 * pow(t, 3);       // N1
            case 2: return h * (t - 2 * pow(t, 2) + pow(t, 3));     // L1
            case 3: return 3 * pow(t, 2) - 2 * pow(t, 3);           // N2
            case 4: return h * (-pow(t, 2) + pow(t, 3));            // L2
            default: break;
          }
        break;

      case 1: // Первые производные
        switch (p)
          {
            case 1: return -6 * t + 6 * pow(t, 2);            // dN1
            case 2: return h * (1 - 4 * t + 3 * pow(t, 2));   // dL1
            case 3: return 6 * t - 6 * pow(t, 2);             // dN2
            case 4: return h * (-2 * t + 3 * pow(t, 2));      // dL2
            default: break;
          }
        break;

      case 2: // Вторые производные
        switch (p)
          {
            case 1: return -6 + 12 * t;                      // ddN1
            case 2: return h * (-4 + 6 * t);                 // ddL1
            case 3: return 6 - 12 * t;                       // ddN2
            case 4: return h * (-2 + 6 * t);                 // ddL2
            default: break;
          }
        break;

      case 3: // Третьи производные
        switch (p)
          {
            case 1: return 12;                               // dddN1
            case 2: return 6 * h;                            // dddL1
            case 3: return -12;                              // dddN2
            case 4: return 6 * h;                            // dddL2
            default: break;
          }
        break;

      default:
        break;
    }

  cout << "OUT OF RANGE base_function: p=" << p << ", derivative=" << derivative << endl;
  return 0.0;
}

/* меняет i, j на номера базисных функций по которым считался интеграл в a_base */
void p_change (int p, int &i, int &j)
{
  switch (p)
    {
      case 0: i = 1; j = 1;   break; 
      case 1: i = 1; j = 2;   break;
      case 2: i = 1; j = 3;   break;
      case 3: i = 1; j = 4;   break;
      case 4: i = 2; j = 2;   break;
      case 5: i = 2; j = 3;   break;
      case 6: i = 2; j = 4;   break;
      case 7: i = 3; j = 3;   break;
      case 8: i = 3; j = 4;   break;
      case 9: i = 4; j = 4;   break;
      default: 
        {
          cout << "OUT OF RANGE p_change" << endl;  break;
        }
    }
}

/* по i,j --> p номер в a_base */
int p_def(int i, int j) 
{
  int p = -1;
  if (i > j)  swap(i, j); // i <= j

  if (i == 1 && j == 1) return 0;  // 11
  if (i == 1 && j == 2) return 1;  // 12
  if (i == 1 && j == 3) return 2;  // 13
  if (i == 1 && j == 4) return 3;  // 14
  if (i == 2 && j == 2) return 4;  // 22
  if (i == 2 && j == 3) return 5;  // 23
  if (i == 2 && j == 4) return 6;  // 24
  if (i == 3 && j == 3) return 7;  // 33
  if (i == 3 && j == 4) return 8;  // 34
  if (i == 4 && j == 4) return 9;  // 44  
  
  return p; 
}

// на каждои отрезке неопределенный инт a_base от 0 до 1 на каждом выбивается h 
double a_base (int p, double h)
{
  switch (p)
    {
      case 0: return  12;              // 36 * t - 72 * pow(t, 2) + 48 * pow(t, 3);             // 11
      case 1: return   6 * h;          // 6 * (4 * t - 7 * pow(t, 2) + 4 * pow(t, 3));          // 12
      case 2: return -12;              // -36 * t + 72 * pow(t, 2) - 48 * pow(t, 3);            // 13
      case 3: return   6 * h;          // 6 * (2 * t - 5 * pow(t, 2) + 4 * pow(t, 3));          // 14
      case 4: return   4 * h * h;      // 16 * t - 24 * pow(t, 2) + 12 * pow(t, 3);             // 22
      case 5: return  -6 * h;          // -6 * (4 * t - 7 * pow(t, 2) + 4 * pow(t, 3));         // 23
      case 6: return   2 * h * h;      // 2 * (4 * t - 9 * pow(t, 2) + 6 * pow(t, 3));          // 24
      case 7: return  12;              // 36 * t - 72 * pow(t, 2) + 48 * pow(t, 3);             // 33
      case 8: return  -6 * h;          // -6 * (2 * t - 5 * pow(t, 2) + 4 * pow(t, 3));         // 34
      case 9: return   4 * h * h;      // 4 * t - 12 * pow(t, 2) + 12 * pow(t, 3);              // 44
      default: 
        {
          cout << "OUT OF RANGE a_base" << endl;   return 0.0;
        }
    }
  return 0;
}

/* считает определенный интеграл от x_i до x_(i+1) и добавляет работу пружины если пружина на K */
double a_res (int p, double x_i, double x_right, double E, double J, double x_k, double k)
{
  double h_cur = x_right - x_i;
  double res = E * J * pow(1 / h_cur, 3) * a_base(p, h_cur);  // 1/h^3 вылезает из вторых производных и якобиана преобразования

  if (x_k < x_i || x_k > x_right)
    {
      return res;
    }
  else // отображаем элемент, на который попал x_k в [0,1] получаем кси_к (определенную от 0 до 1) и на ней считаем произведение 
    {
      int i, j;    p_change(p, i, j);
      return res + k * base_function (i, (x_k - x_i) / h_cur, h_cur, 0) * base_function (j, (x_k - x_i) / h_cur, h_cur, 0);
    }  
}

/* считает определенный интеграл он B_i * phi_j, где B_i линейный базис */
double f_integral (int i, int j, double h)
{
  if (i == 1 && j == 1)   return 3 / 20;
  if (i == 1 && j == 2)   return h * 1 / 30;
  if (i == 1 && j == 3)   return 7 / 20;
  if (i == 1 && j == 4)   return -1 * h / 20;

  if (i == 2 && j == 1)   return 7 / 20;
  if (i == 2 && j == 2)   return h * 1 / 20;
  if (i == 2 && j == 3)   return 3 / 20;
  if (i == 2 && j == 4)   return -1 * h / 30;

  return 0.0;
}

double q_global (double x, double q1, double q2, double x_q1, double x_q2)
{
  if (x < x_q1 || x > x_q2)   return 0.0;

  return q1 + (q2 - q1) / (x_q2 - x_q1) * (x - x_q1);
}
// q2 * t + q1 - q1 * t = (q2 - q1) * t - q1

// подсчет <f,phi_p> с распределенной нагрузкой, моментом и сосредоточенной силой | p - индекс базисных эрмитовых от 1 до 4 
double f_res (int p, double x_i, double x_right, double q1, double q2, double x_q1, double x_q2, double m, double x_m, double p_force, double x_p)
{
  double res = 0.0;
  double h = x_right - x_i;

/* часть с распределенной нагрузкой q 
   по построению сетки конечный элемент может либо лежать целиком на распред. нагрузке, либо целиком не лежать.
   поэтому проверка лежит ли середина правее x_i */
  double median = (x_i + x_right) / 2;
  if (median > x_q1 && median < x_q2)
    {
      //cout << "x_i = " << x_i << "  x_r = " << x_right << "   q_raspred" << endl;
      double q1_, q2_;                              // значения q(x) на концах конечного элемента
      q1_ = q_global(x_i, q1, q2, x_q1, x_q2);
      q2_ = q_global(x_right, q1, q2, x_q1, x_q2);
      res += h * (q2_ * f_integral(1, p, h) + q1_ * f_integral(2, p, h));
    }

/* часть с моментом m */  
  if (x_m >= x_i && x_m < x_right)
    {
      //cout << "x_i = " << x_i << "  x_r = " << x_right << "   moment" << endl;
      res -= m * (1 / h) * base_function (p, (x_m - x_i) / h, h, 1);
    }

/* часть с силой p. Переведенная в локальные координаты */  
  if (x_p >= x_i && x_p < x_right)
    {
      //cout << "x_i = " << x_i << "  x_r = " << x_right << "   force" << endl;
      res += p_force * base_function (p, (x_p - x_i) / h, h, 0);
    }
  return res;  
}

void show_array (double** array, int size)
{
  cout << endl;
  for (int i = 0; i < size; ++i) 
    {
      for (int j = 0; j < size; ++j) 
        {
          cout << setw(6) << setprecision(2) << array[i][j] << " ";
          //cout << setw(12) << array[i][j] << " ";
        }
      cout << endl;
    }
    cout << endl;
}

void array_to_file(const char* filename, double* array, int len)
{
    ofstream outfile(filename);
    if (!outfile.is_open()) {
        cout << "Ошибка: не удалось открыть файл для записи!" << std::endl;
        return;
    }
    int width = 15;
    for (int i = 0; i < len; i++) {
        outfile << scientific << setw (width) << array[i] << endl; 
    }
    outfile.close(); 
}

void show_vect (double* vect, int size, int mode)
{
  for (int i = 0; i < size; ++i) 
    {
      //cout << setw(8) << setprecision(3) << vect[i] << "\n";
      cout  << vect[i] << "  ";
      if (mode == 0)
        {
          cout << endl;
        }
    }
  cout << endl;  
}

void sort (double* array, int size, int* index_v) 
{
  int index_count = 3;
  int original_indices[size];

  for (int i = 0; i < size; i++) 
    original_indices[i] = i;
   
  for (int i = 0; i < size; i++)
  {
    for (int j = i + 1; j < size; j++)
      {
        if (array[i] > array[j])
          {
            swap (array[i], array[j]);
            swap (original_indices[i], original_indices[j]);
          }
      }
  }

  // Обновление массива index_v на основе новых индексов
  for (int i = 0; i < index_count; i++) 
    {
      int original_index = static_cast<int>(index_v[i]);
      
      for (int j = 0; j < size; j++)
        {
          if (original_indices[j] == original_index) 
            {
              index_v[i] = j; // Новый индекс после сортировки
              break;
            }
        }
    }
}

/* 9 обязательных узлов: начало, конец, x_q1, x_r1, x_r2, x_q2, x_r3 x_m x_p / все характерные точки */
int make_mesh (int L, double* data, double* x_mesh, int* index_v) 
{
  if (L < 9)  { cout << "ERROR: L < 9\n";  return -1; }

  /* возвращает h по концам и колву внутреннихъ точек */
  auto into_dots = [](double a, double b, int n) { return fabs(b - a) / static_cast<double>(n + 1); };  

  // [x_k, x_q1, x_q2, x_r1, x_r2, x_r3, x_m, x_z, x_p][q1, q2, p, k, m, E, J]
  x_mesh[L - 9] = 0.0;
  x_mesh[L - 1] = 25.0;

  x_mesh[L - 6] = data[8]; /* x_p */
  x_mesh[L - 2] = data[6]; /* x_m */
  x_mesh[L - 8] = data[1]; /* x_q1 */
  x_mesh[L - 4] = data[2]; /* x_q2 */

  x_mesh[L - 7] = data[3]; /* x_r1 */
  x_mesh[L - 5] = data[4]; /* x_r2 */
  x_mesh[L - 3] = data[5]; /* x_r3 */

  // за этими следим / x_r1  x_r2  x_r3
  index_v[0] = L - 7;    index_v[1] = L - 5;     index_v[2] = L - 3;    index_v[3] = L - 1;

  int min_part = (L - 9) / 8;     // обязательное колво в каждом gate 
  int peace = (L - 9) % 8;        // остаток, раскидываем с хвоста / более равномерно так что ли
  int count = 0;                  // индекс того сколько положили внутренних точек

  for (int i = 0; i < 8; i++)             // i -- обязательныый отрезок
    {
      int n_curr = min_part + ((peace > 0) ? 1 : 0);

      for (size_t j = 1; j <= n_curr; j++)
        {
          x_mesh[count] = x_mesh[L - 2 - i] + into_dots (x_mesh[L - 2 - i], x_mesh[L - 1 - i], n_curr) * j;
          count++;
        }

      if (peace > 0) peace--;
    }

  sort (x_mesh, L, index_v);

  array_to_file ("mesh.txt", x_mesh, L);
  
  return 0;
}

void read_data (const string& filename, double* data)
{
  ifstream inputfile(filename);
  string line;    double number = 0;
  for (int i = 0; i < 16; i++)
    {
      if (getline(inputfile, line))
        {
          istringstream iss(line);
          if (iss >> number)  data[i] = number;
        }
    }
  inputfile.close();
}

// Функция для проверки корректности решения
bool verify_solution (double **A, double *Q, double *F, int size, double tolerance = 1e-9, bool showing_param = 0)
{
  // Вычисляем A * Q
  double *computed_F = new double[size];
  for (int i = 0; i < size; i++)
    {
      computed_F[i] = 0.0;
      for (int j = 0; j < size; j++)
        {
          computed_F[i] += A[i][j] * Q[j];
        }
    }

  for (int i = 0; i < size; i++)  {computed_F[i] = F[i] - computed_F[i]; }
  //cout << "B - A*X = ";   show_vect (computed_F, size, 1);

  // Вычисляем норму разности между F и A * Q
  double norm_diff = 0.0;
  for (int i = 0; i < size; i++)
    {
      norm_diff += computed_F[i] * computed_F[i];
    }
  norm_diff = sqrt(norm_diff);

  // Освобождаем память
  delete[] computed_F;

  // Проверяем, укладывается ли норма в допустимое отклонение
  if (norm_diff < tolerance)
    {
      if (showing_param)
        {
          cout << "Solution verification success: l2_abs norm diff = " << norm_diff << endl;
        }
      return true; // Решение корректно
    }
  else
    {
      if (showing_param)
        {
          cout << "Solution verification failed: l2_abs norm diff = " << norm_diff << endl;
        }
      return false; // Решение некорректно
    }
}

inline void rm_row_column (double** A, double* F, int& size, int rm_index) 
{
  // Удаляем столбец, сдвигая элементы внутри каждой строки
  for (size_t i = 0; i < size; i++) 
    {
      for (size_t j = rm_index; j < size - 1; j++) 
        {
          A[i][j] = A[i][j + 1]; // Сдвиг элементов влево
        }
    }

  // Удаляем строку, сдвигая указатели строк и освобождая лишнюю память
  delete[] A[rm_index]; // Освобождаем память строки, которая больше не используется / потенциально может течь тут

  for (size_t i = rm_index; i < size - 1; i++)    // Сдвиг указателей на строки
    {
      F[i] = F[i + 1];
      A[i] = A[i + 1]; 
    }

  size--;
}

/* index_v = [x_r1, x_r2, x_r3, x_z] */
void rm_index_v (double** A, double* F, int& size, int* index_v)
{
  rm_row_column (A, F, size, 2 * index_v[3] + 1);   // rm last row and column

  for (int i = 3; i > - 1; i--)
    {
      rm_row_column (A, F, size, 2 * index_v[i]);
    }
}

/* добавляет 0 в массив на место index || size = размер укороченного массива */
inline void add_num_to_array (double* array, int& size, int index)
{
  double tmp = array[index];
  array[index] = 0.0;

  for (size_t i = index + 1; i < size + 1; i++)
    {
      swap (tmp, array[i]);
    }

  size++;  
}

inline void add_v_to_solve (double* array, int& size, int* index_v)
{
  for (size_t i = 0; i < 4; i++)
    {
      add_num_to_array (array, size, 2 * index_v[i]);
    }
  add_num_to_array (array, size, 2 * index_v[3] + 1);  
}

void write_arrays_to_file(const char* filename, double* show_mesh, double* show_w, double* show_theta, 
                          double* show_moment, double* show_q, int show_count)
{
    ofstream outfile(filename);
    if (!outfile.is_open()) {
        cout << "Ошибка: не удалось открыть файл для записи!" << std::endl;
        return;
    }

    int width = 15;
    for (int i = 0; i < show_count; i++) {
        outfile << scientific << setw (width) << show_mesh[i] 
                << setw (width) << show_w[i] << " " 
                << setw (width) << show_theta[i] << " " 
                << setw (width) << show_moment[i] << " " 
                << setw (width) << show_q[i] << "\n";
    }

    outfile.close(); 
    //cout << "Данные успешно записаны в файл" << endl;
}

/* строит решения по вектору Q, заполняет массивы w, theta, moment, force */
void make_solution (int show_count, double* show_mesh, double* x_mesh, double* solve, double* w,
                    double* theta, double* moment, double* force, int L, double* data)
{
  double show_h = 25.0 / (show_count - 1);
  for (size_t i = 0; i < show_count; i++)
    {
      show_mesh[i] = i * show_h;
    }

  int mode = 0;
  int count = 0;
  for (int i = 0; i < L - 1; i++)   // по конечным элементам  [x_mesh[i], x_mesh[i + 1]]
    {
      double h_current = x_mesh[i + 1] - x_mesh[i];
      if (h_current < 1e-10)
        cout << "ДЕЛЕНИЕ НА 0 В make_solution" << endl;
      
      if (mode)
        {
          cout << "\ni = " << i << "    h = " << h_current << endl;
        }

      while (show_mesh[count] >= x_mesh[i] && show_mesh[count] <= x_mesh[i + 1])
        {
          if (mode)
            {
              cout << setw (10) << fixed << setprecision(6) << "cur = " << show_mesh[count];
              cout << setw (10) << "       h = " << h_current;
              cout << setw (10) << "       arg = " << (show_mesh[count] - x_mesh[i]) / h_current;
              cout << setw (10) << "       y_1 = " << base_function (1, (show_mesh[count] - x_mesh[i]) / h_current, h_current, 0) << "\n";
            }

          // 3 порядок
          w[count] = solve[2 * i] * base_function (1, (show_mesh[count] - x_mesh[i]) / h_current, h_current, 0) +             
                     solve[2 * i + 1] * base_function (2, (show_mesh[count] - x_mesh[i]) / h_current, h_current, 0) +
                     solve[2 * (i + 1)] * base_function (3, (show_mesh[count] - x_mesh[i]) / h_current, h_current, 0) +
                     solve[2 * (i + 1) + 1] * base_function (4, (show_mesh[count] - x_mesh[i]) / h_current, h_current, 0);
          // 2 порядок
          theta[count] = 1 / h_current * 
                        (solve[2 * i] * base_function (1, (show_mesh[count] - x_mesh[i]) / h_current, h_current, 1) +       
                         solve[2 * i + 1] * base_function (2, (show_mesh[count] - x_mesh[i]) / h_current, h_current, 1) +
                         solve[2 * (i + 1)] * base_function (3, (show_mesh[count] - x_mesh[i]) / h_current, h_current, 1) +
                         solve[2 * (i + 1) + 1] * base_function (4, (show_mesh[count] - x_mesh[i]) / h_current, h_current, 1));
          // 1 порядок
          moment[count] = pow (1 / h_current, 2) * 
                          data[14] * data[15] *
                              (solve[2 * i] * base_function (1, (show_mesh[count] - x_mesh[i]) / h_current, h_current, 2) +
                               solve[2 * i + 1] * base_function (2, (show_mesh[count] - x_mesh[i]) / h_current, h_current, 2) +
                               solve[2 * (i + 1)] * base_function (3, (show_mesh[count] - x_mesh[i]) / h_current, h_current, 2) +
                               solve[2 * (i + 1) + 1] * base_function (4, (show_mesh[count] - x_mesh[i]) / h_current, h_current, 2));   
          // 0 порядок
          force[count] = pow (1 / h_current, 3) *
                        (solve[2 * i] * base_function (1, (show_mesh[count] - x_mesh[i]) / h_current, h_current, 3) +
                         solve[2 * i + 1] * base_function (2, (show_mesh[count] - x_mesh[i]) / h_current, h_current, 3) +
                         solve[2 * (i + 1)] * base_function (3, (show_mesh[count] - x_mesh[i]) / h_current, h_current, 3) +
                         solve[2 * (i + 1) + 1] * base_function (4, (show_mesh[count] - x_mesh[i]) / h_current, h_current, 3)); 

          count++;
        }
    }
}

// Функция для поиска совпадений
void find_indices(double* x_mesh, int mesh_size, double* values, int* found_indices, double tolerance = 1e-6) 
{
  int values_size = 4;

    // Проход по массиву значений
  for (int i = 0; i < values_size; ++i) 
    {
      for (int j = 0; j < mesh_size; ++j) 
        {
          if (abs(x_mesh[j] - values[i]) < tolerance) 
            { 
              found_indices[i] = j;
              break; // Найдено, переходим к следующему значению
          }
      }
  }
}

double* make_uniform_mesh (int &L, int* index_v, int n = 0)
{
  L = 25 * static_cast<int>(pow(2, n)) + 1;
  
  double* x_mesh = new double[L];
  
  double h = 1 / pow (2, n);

  for (int i = 0; i < L; i++)  x_mesh[i] = i * h;

  double values[4] = {4, 12, 16, 25};  
  
  find_indices(x_mesh, L, values, index_v); 
    
  array_to_file("mesh.txt", x_mesh, L);

  return x_mesh;
}