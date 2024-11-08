#include "common.h"
#include "arrays.h"

#include <ccomplex>
#define class_complex std::complex<double>
#define class_create_complex(x,y) std::complex<double>(x,y)
#define creal(x) (x).real()
#define cimag(x) (x).imag()
#define cpow pow
#define cexp exp
#define cabs abs
#define ctan tan
#define csin sin
#define _Complex_I std::complex<double>(0.0,1.0)

#ifdef __cplusplus
extern "C" {
#endif

  int array_integrate_internal_exponential(
                                         double * x0,
                                         int x_size,
                                         int x_stride,
                                         double * y0,
                                         int y_stride,
                                         double * ddy0,
                                         int ddy_stride,
                                         class_complex phase,
                                         class_complex * result,
                                         ErrorMsg errmsg);

  int array_integrate_all_spline_exponential(
            double * array,
            int n_columns,
            int n_lines,
            int index_x,   /** from 0 to (n_columns-1) */
            int index_y,
            int index_ddy,
            class_complex phase,
            class_complex * result,
            ErrorMsg errmsg
            );

  int array_integrate_all_spline_table_lines_exponential(
            double * x,
		        int x_size,
		        double * y_array,
		        int y_size,
		        double * ddy_array,
            class_complex phase,
            class_complex * result,
            ErrorMsg errmsg
            );

  int array_integrate_internal_exponential_pure_phase(
                        double * x0,
                        const int x_size,
                        const int x_stride,
                        double * y0,
                        const int y_stride,
                        double * ddy0,
                        const int ddy_stride,
                        double phase,
                        class_complex * result,
                        class_complex * condition_number,
                        ErrorMsg errmsg);

  int array_integrate_internal_exponential_pure_phase_compensated(
                        double * x0,
                        const int x_size,
                        const int x_stride,
                        double * y0,
                        const int y_stride,
                        double * ddy0,
                        const int ddy_stride,
                        const double phase,
                        class_complex * result,
                        const double condition_number_threshold,
                        class_complex * condition_number,
                        ErrorMsg errmsg);

  int array_integrate_all_spline_fourier_compensated(
            double * array,
            int n_columns,
            int n_lines,
            int index_x,   /** from 0 to (n_columns-1) */
            int index_y,
            int index_ddy,
            double phase,
            class_complex * result,
            const double condition_num_threshold,
            class_complex * condition_num,
            ErrorMsg errmsg
            );

  int array_integrate_all_spline_table_lines_fourier_compensated(
            double * x,
            int x_size,
            double * y_array,
            int y_size,
            double * ddy_array,
            double phase,
            class_complex * result,
            const double condition_num_threshold,
            class_complex * condition_num,
            ErrorMsg errmsg
            );

  int array_square_integrate_all_spline_table_lines(
            double * x,
            int x_size,
            double * y_array,
            int y_size,
            double * ddy_array,
            double * result,
            ErrorMsg errmsg
            );

  int array_square_integrate_exponential_all_spline_table_lines(
            double * x,
            int x_size,
            double * y_array,
            int y_size,
            double * ddy_array,
            const double bias,
            const int derivative_order,
            double * result,
            ErrorMsg errmsg
            );

#ifdef __cplusplus
}
#endif
