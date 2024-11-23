#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <math.h>
#include <cstdio>
#include <cstring>
#include <sstream>
#include <iomanip>
#include <stdexcept> // Для std::runtime_error

using namespace std;

//---------------------------------------------------------------------------------------------------------

void show_vect (double*, int, int);

void show_array (double**, int);

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




/* базисные функции в порядке {N1, L1, N2, L2} */
double base_func (int p, double t)
{
  switch (p)
    {
      case 1: return 1 - 3 * pow(t, 2) + 2 * pow(t, 3);           // N1                           
      case 2: return t - 2 * t * t + pow(t, 3);                   // L1
      case 3: return 3 * pow(t, 2) - 2 * pow(t, 3);               // N2
      case 4: return -t * t + pow(t, 3);                          // L2
      default: 
        {
          cout << "OUT OF RANGE base_func" << endl;   return 0.0;
        }
    }
  return 0;
}

/* производные базисных функций */
double base_der_func (int p, double t)
{
  switch (p)
    {
      case 1: return 6 * t + 6 * pow(t, 2);                   // N1                           
      case 2: return 1 - 4 * t + 3 * pow(t, 2);               // L1
      case 3: return 6 * t - 6 * pow(t, 2);                   // N2
      case 4: return -2 * t + 3 * pow(t, 2);                  // L2
      default: 
        {
          cout << "OUT OF RANGE base_func" << endl;   return 0.0;
        }
    }
  return 0;
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

/* на каждои отрезке неопределенный инт a_base от 0 до 1 */
double a_base (int p, double t)
{
  switch (p)
    {
      case 0: return 36 * t - 72 * pow(t, 2) + 48 * pow(t, 3);                    // 11
      case 1: return 6 * (4 * t - 7 * pow(t, 2) + 4 * pow(t, 3));                 // 12
      case 2: return -36 * t + 72 * pow(t, 2) - 48 * pow(t, 3);                   // 13
      case 3: return 6 * (2 * t - 5 * pow(t, 2) + 4 * pow(t, 3));                 // 14
      case 4: return 16 * t - 24 * pow(t, 2) + 12 * pow(t, 3);                    // 22
      case 5: return -6 * (4 * t - 7 * pow(t, 2) + 4 * pow(t, 3));                // 23
      case 6: return 2 * (4 * t - 9 * pow(t, 2) + 6 * pow(t, 3));                 // 24
      case 7: return 36 * t - 72 * pow(t, 2) + 48 * pow(t, 3);                    // 33
      case 8: return -6 * (2 * t - 5 * pow(t, 2) + 4 * pow(t, 3));                // 34
      case 9: return 4 * t - 12 * pow(t, 2) + 12 * pow(t, 3);                     // 44
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
  double res = E * J * pow(1 / h_cur, 3) * a_base(p, 1.0);  // 1/h^3 вылезает из вторых производных и якобиана преобразования

/* 
  cout << "h_cur = " << h_cur << "  ";
  cout << "p = " << p << "  x_i = " << x_i << "  x_i+1 = " << x_right << "  res = " << res << "   int(1) = " << a_base(p, 1.0);
  cout << "   E*J = " << E*J << "    1/h^3 = " << pow(1 / h_cur, 3) << endl;
*/

  if (x_k < x_i || x_k > x_right)
    {
      return res;
    }
  else /* отображаем элемент, на который попал x_k в [0,1] получаем кси_к (определенную от 0 до 1) и на ней считаем произведение */
    {
      int i,j;    p_change(p, i, j);
      //cout << "  add = " << k * base_func(i, (x_k - x_i) / h_cur) * base_func(j, (x_k - x_i) / h_cur) << endl;
      return res + k * base_func(i, (x_k - x_i) / h_cur) * base_func(j, (x_k - x_i) / h_cur);
    }  
}

/* считает определенный интеграл он B_i * phi_j, где B_i линейный базис */
double f_integral (int i, int j)
{
  if (i == 1 && j == 1)   return 3 / 20;
  if (i == 1 && j == 2)   return 1 / 30;
  if (i == 1 && j == 3)   return 7 / 20;
  if (i == 1 && j == 4)   return -1 / 20;

  if (i == 2 && j == 1)   return 7 / 20;
  if (i == 2 && j == 2)   return 1 / 20;
  if (i == 2 && j == 3)   return 3 / 20;
  if (i == 2 && j == 4)   return -1 / 30;

  return 0.0;
}

double q_global (double x, double q1, double q2, double x_q1, double x_q2)
{
  if (x < x_q1 || x > x_q2)   return 0.0;

  return q1 + (q2 - q1) / (x_q2 - x_q1) * (x - x_q1);
}


/* подсчет <f,phi_p> с распределенной нагрузкой, моментом и сосредоточенной силой 
   p - индекс базисных эрмитовых от 1 до 4 */
double f_res (int p, double x_i, double x_right, double q1, double q2, double x_q1, double x_q2, double m, double x_m, double p_force, double x_p)
{
  double res = 0.0;

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
      res += (x_right - x_i) * (q2_ * f_integral(1, p) + q1_ * f_integral(2, p));
    }


/* часть с моментом m */  
  if (x_m >= x_i && x_m <= x_right)
    {
      //cout << "x_i = " << x_i << "  x_r = " << x_right << "   moment" << endl;
      res += m * base_der_func(p, (x_m - x_i) / (x_right - x_i));
    }

/* часть с силой p. Переведенная в локальные координаты */  
  if (x_p >= x_i && x_p <= x_right)
    {
      //cout << "x_i = " << x_i << "  x_r = " << x_right << "   force" << endl;
      res -= p_force * base_func(p, (x_p - x_i) / (x_right - x_i));
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
          //cout << setw(8) << setprecision(3) << array[i][j] << " ";
          cout << setw(12) << array[i][j] << " ";
        }
      cout << endl;
    }
    cout << endl;
}

void show_vect (double* vect, int size, int mode)
{
  for (int i = 0; i < size; ++i) 
    {
      //cout << setw(8) << setprecision(3) << vect[i] << "\n";
      cout << setw(13) << setprecision(4) << vect[i];
      if (mode == 0)
        {
          cout << endl;
        }
    }
  cout << endl;  
}

void sort (double* array, int size)
{
  for (int i = 0; i < size; i++)
    {
      for (int j = 0; j < size; j++)
        {
          if (array[i] < array[j])
            swap(array[i], array[j]);
        }
    }
}

/* выделяет память под сетку и создает ее */
int make_mesh (int L, double x_q1, double x_q2, double* x_mesh) 
{
  if (L < 4) { cout << "ERROR: L < 4\n";  return -1; }

  double h = 25.0 / double(L - 3);      // равномерный шаг с выбросом точек x_q1 и x_q2

  for (int i = 0; i < L - 2; i++)
    {
      x_mesh[i] = h * i;
    }

  x_mesh[L - 2] = x_q1;    x_mesh[L - 1] = x_q2;
  sort (x_mesh, L);
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
bool verify_solution(double **A, double *Q, double *F, int size, double tolerance = 1e-6)
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

  // Вычисляем норму разности между F и A * Q
  double norm_diff = 0.0;
  for (int i = 0; i < size; i++)
    {
      norm_diff += (F[i] - computed_F[i]) * (F[i] - computed_F[i]);
    }
  norm_diff = sqrt(norm_diff);

  // Освобождаем память
  delete[] computed_F;

  // Проверяем, укладывается ли норма в допустимое отклонение
  if (norm_diff < tolerance)
    {
      return true; // Решение корректно
    }
  else
    {
      cout << "\nSolution verification failed: l2_abs norm difference = " << norm_diff << endl;
      return false; // Решение некорректно
    }
}


// <============================================MAIN==============================================================>
int main (void)
{
  double data[16];    // [x_k, x_q1, x_q2, x_r1, x_r2, x_r3, x_m, x_z, x_p][q1, q2, p, k, m, E, J]
  read_data ("data.txt", data);

//--------------------------------------------------------------------------------- создание сетки   
/* ебло ты утиное */

  int L = 4;            // Матрица A будет размера 2L x 2L.      L - 1 = колво конеч эл.
  double* x_mesh = new double[L];
  if (make_mesh(L, data[1], data[2], x_mesh)) { return -1;} 
  //cout << "x_mesh = ";   show_vect(x_mesh, L, 1);

//---------------------------------------------------------------------------------- выделение памяти под A и F 

  double** A = new double*[2 * L];                                                   
  for (size_t i = 0; i < 2 * L; i++)   
    { 
      A[i] = new double[2 * L];     memset(A[i], 0.0, 2 * L * sizeof(double));
    } 

  double* F = new double[2 * L];    memset(F, 0.0, 2 * L * sizeof(double));  

  double** A_mini = new double*[4];                                                   
  for (size_t i = 0; i < 4; i++)   
    { 
      A_mini[i] = new double[4];
    } 

  double* F_mini = new double[4];   

//-----------------------------------------------------------------------------------

  for (int k = 0; k < L - 1; k++)               // k - индекс конечного элемента
    {
      for (int i = 0; i < 4; i++)
        {
          for (int j = 0; j < 4; j++)
            {
              A[2 * k + i][2 * k + j] += a_res(p_def(i+1, j+1), x_mesh[k], x_mesh[k + 1], data[14], data[15], data[0], data[12]);
              //A_mini[i][j] = a_res(p_def(i+1, j+1), x_mesh[k], x_mesh[k + 1], data[14], data[15], data[0], data[12]);
            }

          F[2 * k + i] += f_res (i + 1, x_mesh[k], x_mesh[k + 1], data[9], data[10], data[1], data[2], data[13], data[6], data[11], data[8]);
          //cout << "k = " << k + 1 << endl;
          //F_mini[i] = f_res (i + 1, x_mesh[k], x_mesh[k + 1], data[9], data[10], data[1], data[2], data[13], data[6], data[11], data[8]);
        }
      //show_array(A_mini, 4);
      //show_vect(F_mini, 4, 0); cout<<endl;
    }


  cout << "\nA и F:\n"; 
  show_array(A, 2 * L);
  show_vect(F, 2 * L, 0);


//-------------------------------------------------------------------------------------


  double *solve = new double[2 * L]; 

  LU_(2 * L, A, F, solve);

  cout << "Q = ";    show_vect(solve, 2 * L, 1);

  verify_solution(A, solve, F, 2 * L);
 
  


//---------------------------------------------------------------------------------------

  delete[] x_mesh;

  return 0;
}