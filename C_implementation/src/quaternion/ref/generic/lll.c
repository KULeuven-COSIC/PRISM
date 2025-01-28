#include <quaternion.h>
#include <rng.h>
#include "internal.h"

// RED(k,l) sub-algorithm
static void
RED(ibz_mat_4x4_t basis, mpf_t u[4][4], mpz_t H[4][4], int k, int l)
{
    mpf_t tmp, tmp2;
    mpz_t q, tmpz;
    mpf_init_set_d(tmp, 0.5);
    mpf_init(tmp2);
    mpz_init(q);
    mpz_init(tmpz);

    // if |u_{k,l}| <= 0.5, terminate
    mpf_abs(tmp2, u[k][l]);
    if (mpf_cmp(tmp2, tmp) <= 0)
        goto end;

    // q <- floor(0.5 + u_{k,l})
    mpf_add(tmp, tmp, u[k][l]);
    mpf_floor(tmp, tmp);
    mpz_set_f(q, tmp);

    // b_k = b_k - q*b_l
    for (int i = 0; i < 4; ++i) {
        mpz_mul(tmpz, q, basis[l][i]);
        mpz_sub(basis[k][i], basis[k][i], tmpz);
    }

    // H_k = H_k - q*H_l
    for (int i = 0; i < 4; ++i) {
        mpz_mul(tmpz, q, H[l][i]);
        mpz_sub(H[k][i], H[k][i], tmpz);
    }

    // u_{k,j} = u_{k,l}-q
    mpf_set_z(tmp2, q);
    mpf_sub(u[k][l], u[k][l], tmp2);

    // forall_i \in 1..l-1: u_{k,i} = u_{k,i} - q*u_{l,i}
    for (int i = 0; i <= l - 1; ++i) {
        mpf_mul(tmp, tmp2, u[l][i]);
        mpf_sub(u[k][i], u[k][i], tmp);
    }

end:
    mpf_clear(tmp);
    mpf_clear(tmp2);
    mpz_clear(q);
    mpz_clear(tmpz);
}

// SWAP(k) sub-algorithm
static void
SWAP(ibz_mat_4x4_t basis,
     mpf_t u[4][4],
     mpz_t H[4][4],
     mpf_t B[4],
     mpf_t bStar[4][4],
     int k,
     int kmax)
{
    mpf_t tmp, tmp2, tmp3, u_tmp, B_tmp, b[4];
    mpf_init(tmp);
    mpf_init(tmp2);
    mpf_init(tmp3);
    mpf_init(u_tmp);
    mpf_init(B_tmp);

    for (int i = 0; i < 4; ++i) {
        mpf_init(b[i]);
    }

    // swap b_k and b_{k-1}
    for (int i = 0; i < 4; ++i) {
        mpz_swap(basis[k][i], basis[k - 1][i]);
    }

    // swap H_k and H_{k-1}
    for (int i = 0; i < 4; ++i) {
        mpz_swap(H[k][i], H[k - 1][i]);
    }

    if (k > 1) {
        // swap u_{k,j} and u_{k-1,j}
        for (int j = 0; j <= k - 2; ++j) {
            mpf_swap(u[k][j], u[k - 1][j]);
        }
    }

    // u = u_{k,k-1}
    mpf_set(u_tmp, u[k][k - 1]);

    // B = B_k + u^2*B_{k-1}
    mpf_mul(B_tmp, u_tmp, u_tmp);
    mpf_mul(B_tmp, B_tmp, B[k - 1]);
    mpf_add(B_tmp, B[k], B_tmp);

    // u_{k,k-1} = u*B_{k-1} / B
    mpf_mul(tmp, u_tmp, B[k - 1]);
    mpf_div(u[k][k - 1], tmp, B_tmp);

    // b = bSTAR_{k-1}
    for (int i = 0; i < 4; ++i) {
        mpf_set(b[i], bStar[k - 1][i]);
    }
    // bSTAR_{k-1}=bSTAR_k+u*b
    for (int i = 0; i < 4; ++i) {
        mpf_mul(tmp, u_tmp, b[i]);
        mpf_add(bStar[k - 1][i], bStar[k][i], tmp);
    }
    // bSTAR_k = -u_{k,k-1}*bSTAR_k+(B_k/B)*b
    mpf_div(tmp2, B[k], B_tmp); // B_k/B
    mpf_neg(tmp, u[k][k - 1]);
    for (int i = 0; i < 4; ++i) {
        mpf_mul(bStar[k][i], tmp, bStar[k][i]);
        mpf_mul(tmp3, tmp2, b[i]);
        mpf_add(bStar[k][i], bStar[k][i], tmp3);
    }

    // B_k = B_{k-1}*B_k/B
    mpf_mul(B[k], B[k - 1], B[k]);
    mpf_div(B[k], B[k], B_tmp);

    // B_{k-1} = B
    mpf_set(B[k - 1], B_tmp);

    for (int i = k + 1; i <= kmax; ++i) {
        // t = u_{i,k}
        mpf_set(tmp, u[i][k]);

        // u_{i,k} = u_{i,k-1} - u*t
        mpf_mul(u[i][k], u_tmp, tmp);
        mpf_sub(u[i][k], u[i][k - 1], u[i][k]);

        // u_{i,k-1} = t + u_{k,k-1}*u_{i,k}
        mpf_mul(tmp2, u[k][k - 1], u[i][k]);
        mpf_add(u[i][k - 1], tmp, tmp2);
    }

    mpf_clear(tmp);
    mpf_clear(tmp2);
    mpf_clear(tmp3);
    mpf_clear(u_tmp);
    mpf_clear(B_tmp);
    for (int i = 0; i < 4; ++i) {
        mpf_clear(b[i]);
    }
}

// m1[0]*m2[0] + m1[1]*m2[1] + q*(m1[2]*m2[2] + m1[3]*m2[3])
static void
dotproduct_row(mpz_t *mul,
               const ibz_mat_4x4_t m1,
               const ibz_mat_4x4_t m2,
               const ibz_t *q,
               int m1j,
               int m2j)
{
    mpz_set_ui(*mul, 0);
    mpz_t tmp1, tmp2;
    mpz_init(tmp1);
    mpz_init(tmp2);
    for (int i = 0; i < 2; ++i) {
        mpz_mul(tmp1, m1[m1j][i], m2[m2j][i]);
        mpz_add(*mul, *mul, tmp1);
    }
    for (int i = 2; i < 4; ++i) {
        mpz_mul(tmp1, m1[m1j][i], m2[m2j][i]);
        mpz_add(tmp2, tmp2, tmp1);
    }
    mpz_mul(tmp2, tmp2, *q);
    mpz_add(*mul, *mul, tmp2);

    mpz_clear(tmp1);
    mpz_clear(tmp2);
}

static void
dotproduct_zr_row(mpf_t *mul,
                  const ibz_mat_4x4_t m1,
                  const mpf_t m2[4][4],
                  const ibz_t *q,
                  int m1j,
                  int m2j)
{
    mpf_set_d(*mul, 0);
    mpf_t tmp1, tmp2;
    mpf_init(tmp1);
    mpf_init(tmp2);
    for (int i = 0; i < 2; ++i) {
        mpf_set_z(tmp1, m1[m1j][i]);
        mpf_mul(tmp1, tmp1, m2[m2j][i]);
        mpf_add(*mul, *mul, tmp1);
    }
    for (int i = 2; i < 4; ++i) {
        mpf_set_z(tmp1, m1[m1j][i]);
        mpf_mul(tmp1, tmp1, m2[m2j][i]);
        mpf_add(tmp2, tmp2, tmp1);
    }
    mpf_set_z(tmp1, *q);
    mpf_mul(tmp2, tmp2, tmp1);
    mpf_add(*mul, *mul, tmp2);

    mpf_clear(tmp1);
    mpf_clear(tmp2);
}

static void
dotproduct_rr_row(mpf_t *mul,
                  const mpf_t m1[4][4],
                  const mpf_t m2[4][4],
                  const ibz_t *q,
                  int m1j,
                  int m2j)
{
    mpf_set_ui(*mul, 0);
    mpf_t tmp1, tmp2;
    mpf_init(tmp1);
    mpf_init(tmp2);
    for (int i = 0; i < 2; ++i) {
        mpf_mul(tmp1, m1[m1j][i], m2[m2j][i]);
        mpf_add(*mul, *mul, tmp1);
    }
    for (int i = 2; i < 4; ++i) {
        mpf_mul(tmp1, m1[m1j][i], m2[m2j][i]);
        mpf_add(tmp2, tmp2, tmp1);
    }
    mpf_set_z(tmp1, *q);
    mpf_mul(tmp2, tmp2, tmp1);
    mpf_add(*mul, *mul, tmp2);

    mpf_clear(tmp1);
    mpf_clear(tmp2);
}

static void
mul_row(mpf_t mul[4][4], const mpf_t *a, const mpf_t m[4][4], int j)
{
    for (int i = 0; i < 4; ++i) {
        mpf_mul(mul[j][i], *a, m[j][i]);
    }
}

static void
add_row(ibz_mat_4x4_t add, const ibz_mat_4x4_t a, const ibz_mat_4x4_t b, int j, int aj, int bj)
{
    for (int i = 0; i < 4; ++i) {
        mpz_add(add[j][i], a[aj][i], b[bj][i]);
    }
}

static void
sub_row(mpf_t add[4][4], const mpf_t a[4][4], const mpf_t b[4][4], int j, int aj, int bj)
{
    for (int i = 0; i < 4; ++i) {
        mpf_sub(add[j][i], a[aj][i], b[bj][i]);
    }
}

/// @brief LLL reduction on 4-dimensional lattice
/// Implements Algorithm 2.6.3 from Henri Cohen's "A Course in Computational Algebraic Number
/// Theory"
/// @param red
/// @param lattice
/// @return
int
quat_lattice_lll(ibz_mat_4x4_t *red, const quat_lattice_t *lattice, const ibz_t *q)
{
    int ret = 0;
    ibz_mat_4x4_t basis;
    mpf_t bStar[4][4];
    mpf_t bStar_tmp[4][4];
    mpf_t tmp;
    mpz_t tmp_z;
    mpf_t cnst;
    mpf_t u[4][4];
    mpz_t H[4][4]; // -> I_4
    mpf_t B[4];

    // Upper bound on log(determinant)
    int logdet = 0;
    for (int i = 0; i < 4; i++) {
      int max = 0;
      for (int j = 0; j < 4; j++) {
	int s = ibz_bitsize(&lattice->basis[i][j]);
	max = s > max ? s : max;
      }
      logdet += max;
    }
    // Set fp precision
    mpf_set_default_prec(2*logdet);
    
    mpf_init(tmp);
    mpz_init(tmp_z);
    mpf_init(cnst);
    for (int i = 0; i < 4; ++i)
        mpf_init(B[i]);

    ibz_mat_4x4_init(&basis);
    ibz_mat_4x4_transpose(&basis, &lattice->basis);

    // Step 1: Initialize: ...
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            mpf_init(u[i][j]);
            mpf_init(bStar[i][j]);
            mpf_init(bStar_tmp[i][j]);
            // bSTAR_1 = b_1 (we copy all)
            if (i == j)
                mpz_init_set_ui(H[i][j], 1);
            else
                mpz_init(H[i][j]);
        }
    }
    int k = 1, kmax = 0;
    // bStar_1 = b_1
    for (int i = 0; i < 4; ++i)
        mpf_set_z(bStar[0][i], basis[0][i]);
    // B_1 = b_1 * b_1
    dotproduct_row(&tmp_z, basis, basis, q, 0, 0);
    mpf_set_z(B[0], tmp_z);

    while (k < 4) {
        // Step 2: Incremental Gram-Schmidt
        // if (k <= kmax) -> we can omit..
        if (k > kmax) {
            kmax = k;
            for (int i = 0; i < 4; ++i) {
                mpf_set_z(bStar[k][i], basis[k][i]);
            }
            for (int j = 0; j <= k - 1; ++j) {
                // bStar_k = b_k -> already done initially -> todo: check if that's ok
                // nop
                // u_{k,j} = b_k*bSTAR_j/B_j
                dotproduct_zr_row(&tmp, basis, bStar, q, k, j);
                mpf_div(u[k][j], tmp, B[j]);
                // bStar_k = bStar_k - u_{k,j}*bStar_j
                mul_row(bStar_tmp, &u[k][j], bStar, j);
                sub_row(bStar, bStar, bStar_tmp, k, k, j);
            }
            // B_k = bStar_k*bStar_k
            dotproduct_rr_row(&B[k], bStar, bStar, q, k, k);
            if (mpf_get_d(B[k]) == 0.0) {
                // b_i did not form a basis, terminate with error
                ret = -1;
                goto err;
            }
        }

        while (1) {
            // Step 3: Test LLL condition
            RED(basis, u, H, k, k - 1);
            // If B_k < (0.75 - u_{k,k-1}^2)*B_{k-1}
            mpf_mul(tmp, u[k][k - 1], u[k][k - 1]);
            mpf_set_d(cnst, 0.99);
            mpf_sub(tmp, cnst, tmp);
            mpf_mul(tmp, tmp, B[k - 1]);
            if (mpf_cmp(B[k], tmp) < 0) {
                SWAP(basis, u, H, B, bStar, k, kmax);
                k = (k - 1 > 1 ? k - 1 : 1);
            } else {
                for (int l = k - 2; l >= 0; --l) {
                    RED(basis, u, H, k, l);
                }
                k++;
                break;
            }
        }
    }
    ibz_mat_4x4_transpose(red, &basis);

err:
    mpf_clear(tmp);
    mpz_clear(tmp_z);
    mpf_clear(cnst);
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            mpf_clear(u[i][j]);
            mpz_clear(H[i][j]);
            mpf_clear(bStar[i][j]);
            mpf_clear(bStar_tmp[i][j]);
        }
    }
    for (int i = 0; i < 4; ++i)
        mpf_clear(B[i]);
    ibz_mat_4x4_finalize(&basis);
    return ret;
}


// #include <quaternion.h>
// #include <rng.h>
// #include "internal.h"

// #if 0

// // RED(k,l) sub-algorithm
// static void RED(ibz_mat_4x4_t basis, mpf_t u[4][4], mpz_t H[4][4], int k, int l) {
//     mpf_t tmp, tmp2;
//     mpz_t q, tmpz;
//     mpf_init_set_d(tmp, 0.5);
//     mpf_init(tmp2);
//     mpz_init(q);
//     mpz_init(tmpz); 

//     // if |u_{k,l}| <= 0.5, terminate
//     mpf_abs(tmp2, u[k][l]);
//     if (mpf_cmp(tmp2, tmp) <= 0)
//         goto end;

//     // q <- floor(0.5 + u_{k,l})
//     mpf_add(tmp, tmp, u[k][l]);
//     mpf_floor(tmp, tmp);
//     mpz_set_f(q, tmp);

//     // b_k = b_k - q*b_l
//     for (int i = 0; i < 4; ++i) {
//         mpz_mul(tmpz, q, basis[l][i]);
//         mpz_sub(basis[k][i], basis[k][i], tmpz);
//     }

//     // H_k = H_k - q*H_l
//     for (int i = 0; i < 4; ++i) {
//         mpz_mul(tmpz, q, H[l][i]);
//         mpz_sub(H[k][i], H[k][i], tmpz);
//     }

//     // u_{k,j} = u_{k,l}-q
//     mpf_set_z(tmp2, q);
//     mpf_sub(u[k][l], u[k][l], tmp2);

//     // forall_i \in 1..l-1: u_{k,i} = u_{k,i} - q*u_{l,i}
//     for (int i = 0; i <= l-1; ++i) {
//         mpf_mul(tmp, tmp2, u[l][i]);
//         mpf_sub(u[k][i], u[k][i], tmp);
//     }

// end:
//     mpf_clear(tmp);
//     mpf_clear(tmp2);
//     mpz_clear(q);
//     mpz_clear(tmpz);
// }

// // SWAP(k) sub-algorithm
// static void SWAP(ibz_mat_4x4_t basis, mpf_t u[4][4], mpz_t H[4][4], mpf_t B[4], mpf_t bStar[4][4], int k, int kmax) {
//     mpf_t tmp, tmp2, tmp3, u_tmp, B_tmp, b[4];
//     mpf_init(tmp);
//     mpf_init(tmp2);
//     mpf_init(tmp3);
//     mpf_init(u_tmp);
//     mpf_init(B_tmp);

//     for (int i = 0; i < 4; ++i) {
//         mpf_init(b[i]);
//     }

//     // swap b_k and b_{k-1}
//     for (int i = 0; i < 4; ++i) {
//         mpz_swap(basis[k][i], basis[k-1][i]);
//     }

//     // swap H_k and H_{k-1}
//     for (int i = 0; i < 4; ++i) {
//         mpz_swap(H[k][i], H[k-1][i]);
//     }

//     if (k > 1) {
//         // swap u_{k,j} and u_{k-1,j}
//         for (int j = 0; j <= k - 2; ++j) {
//             mpf_swap(u[k][j], u[k-1][j]);
//         }
//     }

//     // u = u_{k,k-1}
//     mpf_set(u_tmp, u[k][k - 1]);

//     // B = B_k + u^2*B_{k-1}
//     mpf_mul(B_tmp, u_tmp, u_tmp);
//     mpf_mul(B_tmp, B_tmp, B[k-1]);
//     mpf_add(B_tmp, B[k], B_tmp);

//     // u_{k,k-1} = u*B_{k-1} / B
//     mpf_mul(tmp, u_tmp, B[k-1]);
//     mpf_div(u[k][k-1], tmp, B_tmp);

//     // b = bSTAR_{k-1}
//     for (int i = 0; i < 4; ++i) {
//         mpf_set(b[i], bStar[k-1][i]);
//     }
//     // bSTAR_{k-1}=bSTAR_k+u*b
//     for (int i = 0; i < 4; ++i) {
//         mpf_mul(tmp, u_tmp, b[i]);
//         mpf_add(bStar[k-1][i], bStar[k][i], tmp);
//     }
//     // bSTAR_k = -u_{k,k-1}*bSTAR_k+(B_k/B)*b
//     mpf_div(tmp2, B[k], B_tmp); // B_k/B
//     mpf_neg(tmp, u[k][k-1]);
//     for (int i = 0; i < 4; ++i) {
//         mpf_mul(bStar[k][i], tmp, bStar[k][i]);
//         mpf_mul(tmp3, tmp2, b[i]);
//         mpf_add(bStar[k][i], bStar[k][i], tmp3);
//     }

//     // B_k = B_{k-1}*B_k/B
//     mpf_mul(B[k], B[k-1], B[k]);
//     mpf_div(B[k], B[k], B_tmp);

//     // B_{k-1} = B
//     mpf_set(B[k-1], B_tmp);

//     for (int i = k+1; i <= kmax; ++i) {
//         // t = u_{i,k}
//         mpf_set(tmp, u[i][k]);

//         // u_{i,k} = u_{i,k-1} - u*t
//         mpf_mul(u[i][k], u_tmp, tmp);
//         mpf_sub(u[i][k], u[i][k-1], u[i][k]);

//         // u_{i,k-1} = t + u_{k,k-1}*u_{i,k}
//         mpf_mul(tmp2, u[k][k-1], u[i][k]);
//         mpf_add(u[i][k-1], tmp, tmp2);
//     }

//     mpf_clear(tmp);
//     mpf_clear(tmp2);
//     mpf_clear(tmp3);
//     mpf_clear(u_tmp);
//     mpf_clear(B_tmp);
//     for (int i = 0; i < 4; ++i) {
//         mpf_clear(b[i]);
//     }

// }

// // m1[0]*m2[0] + m1[1]*m2[1] + q*(m1[2]*m2[2] + m1[3]*m2[3])
// static void dotproduct_row(mpz_t* mul, const ibz_mat_4x4_t m1, const ibz_mat_4x4_t m2, const ibz_t *q, int m1j, int m2j) {
//     mpz_set_ui(*mul, 0);
//     mpz_t tmp1, tmp2;
//     mpz_init(tmp1);
//     mpz_init(tmp2);
//     for (int i = 0; i < 2; ++i) {
//         mpz_mul(tmp1, m1[m1j][i], m2[m2j][i]);
//         mpz_add(*mul, *mul, tmp1);
//     }
//     for (int i = 2; i < 4; ++i) {
//         mpz_mul(tmp1, m1[m1j][i], m2[m2j][i]);
//         mpz_add(tmp2, tmp2, tmp1);
//     }
//     mpz_mul(tmp2, tmp2, *q);
//     mpz_add(*mul, *mul, tmp2);

//     mpz_clear(tmp1);
//     mpz_clear(tmp2);
// }

// static void dotproduct_zr_row(mpf_t* mul, const ibz_mat_4x4_t m1, const mpf_t m2[4][4], const ibz_t *q, int m1j, int m2j) {
//     mpf_set_d(*mul, 0);
//     mpf_t tmp1, tmp2;
//     mpf_init(tmp1);
//     mpf_init(tmp2);
//     for (int i = 0; i < 2; ++i) {
//         mpf_set_z(tmp1, m1[m1j][i]);
//         mpf_mul(tmp1, tmp1, m2[m2j][i]);
//         mpf_add(*mul, *mul, tmp1); 
//     }
//     for (int i = 2; i < 4; ++i) {
//         mpf_set_z(tmp1, m1[m1j][i]);
//         mpf_mul(tmp1, tmp1, m2[m2j][i]);
//         mpf_add(tmp2, tmp2, tmp1); 
//     }
//     mpf_set_z(tmp1, *q);
//     mpf_mul(tmp2, tmp2, tmp1);
//     mpf_add(*mul, *mul, tmp2);

//     mpf_clear(tmp1);
//     mpf_clear(tmp2);
// }

// static void dotproduct_rr_row(mpf_t* mul, const mpf_t m1[4][4], const mpf_t m2[4][4], const ibz_t *q, int m1j, int m2j) {
//     mpf_set_ui(*mul, 0);
//     mpf_t tmp1, tmp2;
//     mpf_init(tmp1);
//     mpf_init(tmp2);
//     for (int i = 0; i < 2; ++i) {
//         mpf_mul(tmp1, m1[m1j][i], m2[m2j][i]);
//         mpf_add(*mul, *mul, tmp1);
//     }
//     for (int i = 2; i < 4; ++i) {
//         mpf_mul(tmp1, m1[m1j][i], m2[m2j][i]);
//         mpf_add(tmp2, tmp2, tmp1);
//     }
//     mpf_set_z(tmp1, *q);
//     mpf_mul(tmp2, tmp2, tmp1);
//     mpf_add(*mul, *mul, tmp2);

//     mpf_clear(tmp1);
//     mpf_clear(tmp2);
// }

// static void mul_row(mpf_t mul[4][4], const mpf_t* a, const mpf_t m[4][4], int j) {
//     for (int i = 0; i < 4; ++i) {
//         mpf_mul(mul[j][i], *a, m[j][i]);
//     }
// }

// static void add_row(ibz_mat_4x4_t add, const ibz_mat_4x4_t a, const ibz_mat_4x4_t b, int j, int aj, int bj) {
//     for (int i = 0; i < 4; ++i) {
//         mpz_add(add[j][i], a[aj][i], b[bj][i]);
//     }
// }

// static void sub_row(mpf_t add[4][4], const mpf_t a[4][4], const mpf_t b[4][4], int j, int aj, int bj) {
//     for (int i = 0; i < 4; ++i) {
//         mpf_sub(add[j][i], a[aj][i], b[bj][i]);
//     }
// }

// /// @brief LLL reduction on 4-dimensional lattice
// /// Implements Algorithm 2.6.3 from Henri Cohen's "A Course in Computational Algebraic Number Theory"
// /// @param red 
// /// @param lattice 
// /// @return 
// int quat_lattice_lll(ibz_mat_4x4_t *red, const quat_lattice_t *lattice, const ibz_t *q, int precision) {
//     if (precision != 0)
//         mpf_set_default_prec(precision);
//     int ret = 0;
//     ibz_mat_4x4_t basis;
//     mpf_t bStar[4][4];
//     mpf_t bStar_tmp[4][4];
//     mpf_t tmp;
//     mpz_t tmp_z;
//     mpf_t cnst;
//     mpf_t u[4][4];
//     mpz_t H[4][4]; // -> I_4
//     mpf_t B[4];
//     mpf_init(tmp);
//     mpz_init(tmp_z);
//     mpf_init(cnst);
//     for (int i = 0; i < 4; ++i)
//         mpf_init(B[i]);

//     ibz_mat_4x4_init(&basis);
//     ibz_mat_4x4_transpose(&basis, &lattice->basis);

//     // Step 1: Initialize: ...
//     for (int i = 0; i < 4; ++i) {
//         for (int j = 0; j < 4; ++j) {
//             mpf_init(u[i][j]);
//             mpf_init(bStar[i][j]);
//             mpf_init(bStar_tmp[i][j]);
//             // bSTAR_1 = b_1 (we copy all)
//             if (i == j)
//                 mpz_init_set_ui(H[i][j], 1);
//             else
//                 mpz_init(H[i][j]);
//         }
//     }
//     int k = 1, kmax = 0;
//     // bStar_1 = b_1
//     for (int i = 0; i < 4; ++i)
//         mpf_set_z(bStar[0][i], basis[0][i]);
//     // B_1 = b_1 * b_1
//     dotproduct_row(&tmp_z, basis, basis, q, 0, 0);
//     mpf_set_z(B[0], tmp_z);

//     while (k < 4) {
//         // Step 2: Incremental Gram-Schmidt
//         // if (k <= kmax) -> we can omit..
//         if (k > kmax) {
//             kmax = k;
//             for (int i = 0; i < 4; ++i) {
//                 mpf_set_z(bStar[k][i], basis[k][i]);
//             }
//             for (int j = 0; j <= k-1; ++j) {
//                 // bStar_k = b_k -> already done initially -> todo: check if that's ok
//                 // nop
//                 // u_{k,j} = b_k*bSTAR_j/B_j
//                 dotproduct_zr_row(&tmp, basis, bStar, q, k, j);
//                 mpf_div(u[k][j], tmp, B[j]);
//                 // bStar_k = bStar_k - u_{k,j}*bStar_j
//                 mul_row(bStar_tmp, &u[k][j], bStar, j);
//                 sub_row(bStar, bStar, bStar_tmp, k, k, j);
//             }
//             // B_k = bStar_k*bStar_k
//             dotproduct_rr_row(&B[k], bStar, bStar, q, k, k);
//             if (mpf_get_d(B[k]) == 0.0) {
//                 // b_i did not form a basis, terminate with error
//                 ret = -1;
//                 goto err;
//             }
//         }

//         while(1) {
//             // Step 3: Test LLL condition
//             RED(basis, u, H, k, k - 1);
//             // If B_k < (0.75 - u_{k,k-1}^2)*B_{k-1}
//             mpf_mul(tmp, u[k][k-1], u[k][k-1]);
//             mpf_set_d(cnst, 0.99);
//             mpf_sub(tmp, cnst, tmp);
//             mpf_mul(tmp, tmp, B[k-1]);
//             if (mpf_cmp(B[k], tmp) < 0) {
//                 SWAP(basis, u, H, B, bStar, k, kmax);
//                 k = (k - 1 > 1 ? k - 1 : 1);
//             } else {
//                 for (int l = k - 2; l >= 0; --l) {
//                     RED(basis, u, H, k, l);
//                 }
//                 k++;
//                 break;
//             }
//         }
//     }
//     ibz_mat_4x4_transpose(red, &basis);

// err:
//     mpf_clear(tmp);
//     mpz_clear(tmp_z);
//     mpf_clear(cnst);
//     for (int i = 0; i < 4; ++i) {
//         for (int j = 0; j < 4; ++j) {
//             mpf_clear(u[i][j]);
//             mpz_clear(H[i][j]);
//             mpf_clear(bStar[i][j]);
//             mpf_clear(bStar_tmp[i][j]);
//         }
//     }
//     for (int i = 0; i < 4; ++i)
//         mpf_clear(B[i]);
//     ibz_mat_4x4_finalize(&basis);
//     return ret;
// }

// #else 


// // Changes compared to floating-point:
// // ibz_div_floor can be moved to HNF
// // rationals only need to be compiled for tests
// // no need to include interals
// // RED(k,l) sub-algorithm
// static void
// RED(ibz_mat_4x4_t *basis, ibq_t (*u)[4][4], ibz_t (*H)[4][4], int k, int l)
// {
//     ibq_t tmp, tmp2;
//     ibz_t q, tibz, num, den, r;
//     ibq_init(&tmp);
//     ibq_init(&tmp2);
//     ibz_init(&q);
//     ibz_init(&tibz);
//     ibz_init(&num);
//     ibz_init(&den);
//     ibz_init(&r);
//     ibz_set(&num, 1);
//     ibz_set(&den, 2);
//     ibq_set(&tmp, &num, &den);
//     ibz_set(&num, 0);
//     ibz_set(&den, 0);

//     // if |u_{k,l}| <= 0.5, terminate
//     ibq_abs(&tmp2, &((*u)[k][l]));
//     if (ibq_cmp(&tmp2, &tmp) <= 0)
//         goto end;

//     // q <- floor(0.5 + u_{k,l})
//     ibq_add(&tmp, &tmp, &((*u)[k][l]));

//     ibq_num(&num, &tmp);
//     ibq_denom(&den, &tmp);
//     // FDIV was used, needs reeimplementation
//     ibz_div_floor(&q, &r, &num, &den);
//     // ibq_floor(tmp, tmp);
//     // ibz_set_f(q, tmp);

//     // b_k = b_k - q*b_l
//     for (int i = 0; i < 4; ++i) {
//         ibz_mul(&tibz, &q, &((*basis)[l][i]));
//         ibz_sub(&((*basis)[k][i]), &((*basis)[k][i]), &tibz);
//     }

//     // H_k = H_k - q*H_l
//     for (int i = 0; i < 4; ++i) {
//         ibz_mul(&tibz, &q, &((*H)[l][i]));
//         ibz_sub(&((*H)[k][i]), &((*H)[k][i]), &tibz);
//     }

//     // u_{k,j} = u_{k,l}-q
//     ibq_set(&tmp2, &q, &ibz_const_one);
//     ibq_sub(&((*u)[k][l]), &((*u)[k][l]), &tmp2);

//     // forall_i \in 1..l-1: u_{k,i} = u_{k,i} - q*u_{l,i}
//     for (int i = 0; i <= l - 1; ++i) {
//         ibq_mul(&tmp, &tmp2, &((*u)[l][i]));
//         ibq_sub(&((*u)[k][i]), &((*u)[k][i]), &tmp);
//     }

// end:
//     ibq_finalize(&tmp);
//     ibq_finalize(&tmp2);
//     ibz_finalize(&q);
//     ibz_finalize(&tibz);
//     ibz_finalize(&num);
//     ibz_finalize(&den);
//     ibz_finalize(&r);
// }

// // SWAP(k) sub-algorithm
// static void
// SWAP(ibz_mat_4x4_t *basis,
//      ibq_t (*u)[4][4],
//      ibz_t (*H)[4][4],
//      ibq_t (*B)[4],
//      ibq_t (*bStar)[4][4],
//      int k,
//      int kmax)
// {
//     ibq_t tmp, tmp2, tmp3, u_tmp, B_tmp, b[4];
//     ibq_init(&tmp);
//     ibq_init(&tmp2);
//     ibq_init(&tmp3);
//     ibq_init(&u_tmp);
//     ibq_init(&B_tmp);

//     for (int i = 0; i < 4; ++i) {
//         ibq_init(&(b[i]));
//     }

//     // swap b_k and b_{k-1}
//     for (int i = 0; i < 4; ++i) {
//         ibz_swap(&((*basis)[k][i]), &((*basis)[k - 1][i]));
//     }

//     // swap H_k and H_{k-1}
//     for (int i = 0; i < 4; ++i) {
//         ibz_swap(&((*H)[k][i]), &((*H)[k - 1][i]));
//     }

//     if (k > 1) {
//         // swap u_{k,j} and u_{k-1,j}
//         for (int j = 0; j <= k - 2; ++j) {
//             ibq_swap(&((*u)[k][j]), &((*u)[k - 1][j]));
//         }
//     }

//     // u = u_{k,k-1}
//     ibq_copy(&u_tmp, &((*u)[k][k - 1]));

//     // B = B_k + u^2*B_{k-1}
//     ibq_mul(&B_tmp, &u_tmp, &u_tmp);
//     ibq_mul(&B_tmp, &B_tmp, &((*B)[k - 1]));
//     ibq_add(&B_tmp, &((*B)[k]), &B_tmp);

//     // u_{k,k-1} = u*B_{k-1} / B
//     ibq_mul(&tmp, &u_tmp, &((*B)[k - 1]));
//     ibq_div(&((*u)[k][k - 1]), &tmp, &B_tmp);

//     // b = bSTAR_{k-1}
//     for (int i = 0; i < 4; ++i) {
//         ibq_copy(&(b[i]), &((*bStar)[k - 1][i]));
//     }
//     // bSTAR_{k-1}=bSTAR_k+u*b
//     for (int i = 0; i < 4; ++i) {
//         ibq_mul(&tmp, &u_tmp, &(b[i]));
//         ibq_add(&((*bStar)[k - 1][i]), &((*bStar)[k][i]), &tmp);
//     }
//     // bSTAR_k = -u_{k,k-1}*bSTAR_k+(B_k/B)*b
//     ibq_div(&tmp2, &((*B)[k]), &B_tmp); // B_k/B
//     ibq_neg(&tmp, &((*u)[k][k - 1]));
//     for (int i = 0; i < 4; ++i) {
//         ibq_mul(&((*bStar)[k][i]), &tmp, &((*bStar)[k][i]));
//         ibq_mul(&tmp3, &tmp2, &(b[i]));
//         ibq_add(&((*bStar)[k][i]), &((*bStar)[k][i]), &tmp3);
//     }

//     // B_k = B_{k-1}*B_k/B
//     ibq_mul(&((*B)[k]), &((*B)[k - 1]), &((*B)[k]));
//     ibq_div(&((*B)[k]), &((*B)[k]), &B_tmp);

//     // B_{k-1} = B
//     ibq_copy(&((*B)[k - 1]), &B_tmp);

//     for (int i = k + 1; i <= kmax; ++i) {
//         // t = u_{i,k}
//         ibq_copy(&tmp, &((*u)[i][k]));

//         // u_{i,k} = u_{i,k-1} - u*t
//         ibq_mul(&((*u)[i][k]), &u_tmp, &tmp);
//         ibq_sub(&((*u)[i][k]), &((*u)[i][k - 1]), &((*u)[i][k]));

//         // u_{i,k-1} = t + u_{k,k-1}*u_{i,k}
//         ibq_mul(&tmp2, &((*u)[k][k - 1]), &((*u)[i][k]));
//         ibq_add(&((*u)[i][k - 1]), &tmp, &tmp2);
//     }

//     ibq_finalize(&tmp);
//     ibq_finalize(&tmp2);
//     ibq_finalize(&tmp3);
//     ibq_finalize(&u_tmp);
//     ibq_finalize(&B_tmp);
//     for (int i = 0; i < 4; ++i) {
//         ibq_finalize(&(b[i]));
//     }
// }

// // m1[0]*m2[0] + m1[1]*m2[1] + q*(m1[2]*m2[2] + m1[3]*m2[3])
// static void
// dotproduct_row(ibz_t *mul,
//                const ibz_mat_4x4_t *m1,
//                const ibz_mat_4x4_t *m2,
//                const ibz_t *q,
//                int m1j,
//                int m2j)
// {
//     ibz_set(mul, 0);
//     ibz_t tmp1, tmp2;
//     ibz_init(&tmp1);
//     ibz_init(&tmp2);
//     for (int i = 0; i < 2; ++i) {
//         ibz_mul(&tmp1, &((*m1)[m1j][i]), &((*m2)[m2j][i]));
//         ibz_add(mul, mul, &tmp1);
//     }
//     for (int i = 2; i < 4; ++i) {
//         ibz_mul(&tmp1, &((*m1)[m1j][i]), &((*m2)[m2j][i]));
//         ibz_add(&tmp2, &tmp2, &tmp1);
//     }
//     ibz_mul(&tmp2, &tmp2, q);
//     ibz_add(mul, mul, &tmp2);

//     ibz_finalize(&tmp1);
//     ibz_finalize(&tmp2);
// }

// static void
// dotproduct_zr_row(ibq_t *mul,
//                   const ibz_mat_4x4_t *m1,
//                   const ibq_t (*m2)[4][4],
//                   const ibz_t *q,
//                   int m1j,
//                   int m2j)
// {
//     ibq_set(mul, &ibz_const_zero, &ibz_const_one);
//     ibq_t tmp1, tmp2;
//     ibq_init(&tmp1);
//     ibq_init(&tmp2);
//     for (int i = 0; i < 2; ++i) {
//         ibq_set(&tmp1, &((*m1)[m1j][i]), &ibz_const_one);
//         ibq_mul(&tmp1, &tmp1, &((*m2)[m2j][i]));
//         ibq_add(mul, mul, &tmp1);
//     }
//     for (int i = 2; i < 4; ++i) {
//         ibq_set(&tmp1, &((*m1)[m1j][i]), &ibz_const_one);
//         ibq_mul(&tmp1, &tmp1, &((*m2)[m2j][i]));
//         ibq_add(&tmp2, &tmp2, &tmp1);
//     }
//     ibq_set(&tmp1, q, &ibz_const_one);
//     ibq_mul(&tmp2, &tmp2, &tmp1);
//     ibq_add(mul, mul, &tmp2);

//     ibq_finalize(&tmp1);
//     ibq_finalize(&tmp2);
// }

// static void
// dotproduct_rr_row(ibq_t *mul,
//                   const ibq_t (*m1)[4][4],
//                   const ibq_t (*m2)[4][4],
//                   const ibz_t *q,
//                   int m1j,
//                   int m2j)
// {
//     // ibq_set(mul, 0);
//     ibq_set(mul, &ibz_const_zero, &ibz_const_one);
//     ibq_t tmp1, tmp2;
//     ibq_init(&tmp1);
//     ibq_init(&tmp2);
//     for (int i = 0; i < 2; ++i) {
//         ibq_mul(&tmp1, &((*m1)[m1j][i]), &((*m2)[m2j][i]));
//         ibq_add(mul, mul, &tmp1);
//     }
//     for (int i = 2; i < 4; ++i) {
//         ibq_mul(&tmp1, &((*m1)[m1j][i]), &((*m2)[m2j][i]));
//         ibq_add(&tmp2, &tmp2, &tmp1);
//     }
//     ibq_set(&tmp1, q, &ibz_const_one);
//     ibq_mul(&tmp2, &tmp2, &tmp1);
//     ibq_add(mul, mul, &tmp2);

//     ibq_finalize(&tmp1);
//     ibq_finalize(&tmp2);
// }

// static void
// mul_row(ibq_t (*mul)[4][4], const ibq_t *a, const ibq_t (*m)[4][4], int j)
// {
//     for (int i = 0; i < 4; ++i) {
//         ibq_mul(&((*mul)[j][i]), a, &((*m)[j][i]));
//     }
// }

// static void
// add_row(ibz_mat_4x4_t *add, const ibz_mat_4x4_t *a, const ibz_mat_4x4_t *b, int j, int aj, int bj)
// {
//     for (int i = 0; i < 4; ++i) {
//         ibz_add(&((*add)[j][i]), &((*a)[aj][i]), &((*b)[bj][i]));
//     }
// }

// static void
// sub_row(ibq_t (*add)[4][4], const ibq_t (*a)[4][4], const ibq_t (*b)[4][4], int j, int aj, int bj)
// {
//     for (int i = 0; i < 4; ++i) {
//         ibq_sub(&((*add)[j][i]), &((*a)[aj][i]), &((*b)[bj][i]));
//     }
// }

// /// @brief LLL reduction on 4-dimensional lattice
// /// Implements Algorithm 2.6.3 from Henri Cohen's "A Course in Computational Algebraic Number
// /// Theory"
// /// @param red
// /// @param lattice
// /// @return
// int
// quat_lattice_lll(ibz_mat_4x4_t *red, const quat_lattice_t *lattice, const ibz_t *p, int precision)
// {
//     if (precision != 0)
//         mpf_set_default_prec(precision);
//     int ret = 0;
//     ibz_mat_4x4_t basis;
//     ibq_t bStar[4][4];
//     ibq_t bStar_tmp[4][4];
//     ibq_t tmp;
//     ibz_t tmp_z;
//     ibz_t den;
//     ibz_t num;
//     ibq_t cnst;
//     ibq_t u[4][4];
//     ibz_t H[4][4]; // -> I_4
//     ibq_t B[4];
//     ibq_init(&tmp);
//     ibz_init(&tmp_z);
//     ibz_init(&den);
//     ibz_init(&num);
//     ibq_init(&cnst);
//     for (int i = 0; i < 4; ++i)
//         ibq_init(&(B[i]));

//     ibz_mat_4x4_init(&basis);
//     ibz_mat_4x4_transpose(&basis, &(lattice->basis));

//     // Step 1: Initialize: ...
//     for (int i = 0; i < 4; ++i) {
//         for (int j = 0; j < 4; ++j) {
//             ibq_init(&(u[i][j]));
//             ibq_init(&(bStar[i][j]));
//             ibq_init(&(bStar_tmp[i][j]));
//             // bSTAR_1 = b_1 (we copy all)
//             if (i == j) {
//                 ibz_init(&(H[i][j]));
//                 ibz_set(&(H[i][j]), 1);
//             } else {
//                 ibz_init(&(H[i][j]));
//             }
//         }
//     }
//     int k = 1, kmax = 0;
//     // bStar_1 = b_1
//     for (int i = 0; i < 4; ++i)
//         ibq_set(&(bStar[0][i]), &(basis[0][i]), &ibz_const_one);
//     // B_1 = b_1 * b_1
//     dotproduct_row(&tmp_z, &basis, &basis, p, 0, 0);
//     ibq_set(&(B[0]), &tmp_z, &ibz_const_one);
//     ibz_set(&num, 99);
//     ibz_set(&den, 100);
//     ibq_set(&cnst, &num, &den);

//     while (k < 4) {
//         // Step 2: Incremental Gram-Schmidt
//         // if (k <= kmax) -> we can omit..
//         if (k > kmax) {
//             kmax = k;
//             for (int i = 0; i < 4; ++i) {
//                 ibq_set(&(bStar[k][i]), &(basis[k][i]), &ibz_const_one);
//             }
//             for (int j = 0; j <= k - 1; ++j) {
//                 // bStar_k = b_k -> already done initially
//                 // nop
//                 // u_{k,j} = b_k*bSTAR_j/B_j
//                 dotproduct_zr_row(&tmp, &basis, &bStar, p, k, j);
//                 ibq_div(&(u[k][j]), &tmp, &(B[j]));
//                 // bStar_k = bStar_k - u_{k,j}*bStar_j
//                 mul_row(&bStar_tmp, &(u[k][j]), &bStar, j);
//                 sub_row(&bStar, &bStar, &bStar_tmp, k, k, j);
//             }
//             // B_k = bStar_k*bStar_k
//             dotproduct_rr_row(&(B[k]), &bStar, &bStar, p, k, k);
//             if (ibq_is_zero(&(B[k]))) {
//                 // b_i did not form a basis, terminate with error
//                 ret = -1;
//                 goto err;
//             }
//         }

//         while (1) {
//             // Step 3: Test LLL condition
//             RED(&basis, &u, &H, k, k - 1);
//             // If B_k < (0.75 - u_{k,k-1}^2)*B_{k-1}
//             ibq_mul(&tmp, &(u[k][k - 1]), &(u[k][k - 1]));
//             ibq_sub(&tmp, &cnst, &tmp);
//             ibq_mul(&tmp, &tmp, &(B[k - 1]));
//             if (ibq_cmp(&(B[k]), &tmp) < 0) {
//                 SWAP(&basis, &u, &H, &B, &bStar, k, kmax);
//                 k = (k - 1 > 1 ? k - 1 : 1);
//             } else {
//                 for (int l = k - 2; l >= 0; --l) {
//                     RED(&basis, &u, &H, k, l);
//                 }
//                 k++;
//                 break;
//             }
//         }
//     }
//     ibz_mat_4x4_transpose(red, &basis);

// err:
//     ibq_finalize(&tmp);
//     ibz_finalize(&tmp_z);
//     ibz_finalize(&num);
//     ibz_finalize(&den);
//     ibq_finalize(&cnst);
//     for (int i = 0; i < 4; ++i) {
//         for (int j = 0; j < 4; ++j) {
//             ibq_finalize(&(u[i][j]));
//             ibz_finalize(&(H[i][j]));
//             ibq_finalize(&(bStar[i][j]));
//             ibq_finalize(&(bStar_tmp[i][j]));
//         }
//     }
//     for (int i = 0; i < 4; ++i)
//         ibq_finalize(&(B[i]));
//     ibz_mat_4x4_finalize(&basis);
//     return ret;
// }



// #endif 




// #if 0
// // double versions

// // RED(k,l) sub-algorithm
// static void RED_dbl(int basis[4][4], double u[4][4], int H[4][4], int k, int l) {

//     // if |u_{k,l}| <= 0.5, terminate
//     if (fabs(u[k][l]) <= 0.5)
//         return;

//     // q <- floor(0.5 + u_{k,l})
//     int q = (int)floor(0.5 + u[k][l]);

//     // b_k = b_k - q*b_l
//     for (int i = 0; i < 4; ++i) {
//         basis[k][i] = basis[k][i] - q*basis[l][i];
//     }

//     // H_k = H_k - q*H_l
//     for (int i = 0; i < 4; ++i) {
//         H[k][i] = H[k][i] - q*H[l][i];
//     }

//     // u_{k,l} = u_{k,l}-q
//     u[k][l] = u[k][l] - q;

//     // forall_i \in 1..l-1: u_{k,i} = u_{k,i} - q*u_{l,i}
//     for (int i = 0; i <= l-1; ++i) { // check: of i < l - 1
//         u[k][i] = u[k][i] - q*u[l][i];
//     }

// }

// // SWAP(k) sub-algorithm
// static void SWAP_dbl(int basis[4][4], double u[4][4], int H[4][4], double B[4], double bStar[4][4], int k, int kmax) {
//     double b[4] = {0};

//     // swap b_k and b_{k-1}
//     for (int i = 0; i < 4; ++i) {
//         int tmp = basis[k][i];
//         basis[k][i] = basis[k-1][i];
//         basis[k-1][i] = tmp;
//     }

//     // swap H_k and H_{k-1}
//     for (int i = 0; i < 4; ++i) {
//         int tmp = H[k][i];
//         H[k][i] = H[k-1][i];
//         H[k-1][i] = tmp;
//     }

//     if (k > 1) {
//         // swap u_{k,j} and u_{k-1,j}
//         for (int j = 0; j <= k - 2; ++j) {
//             double tmp = u[k][j];
//             u[k][j] = u[k-1][j];
//             u[k-1][j] = tmp;
//         }
//     }

//     // u = u_{k,k-1}
//     double u_tmp = u[k][k-1];

//     // B = B_k + u^2*B_{k-1}
//     double B_tmp = B[k] + u_tmp*u_tmp*B[k-1];

//     // u_{k,k-1} = u*B_{k-1} / B
//     u[k][k-1] = u_tmp * B[k-1] / B_tmp;

//     // b = bSTAR_{k-1}
//     for (int i = 0; i < 4; ++i) {
//         b[i] = bStar[k-1][i];
//     }
//     // bSTAR_{k-1}=bSTAR_k+u*b
//     for (int i = 0; i < 4; ++i) {
//         bStar[k-1][i] = bStar[k][i]+u_tmp*b[i];
//     }
//     // bSTAR_k = -u_{k,k-1}*bSTAR_k+(B_k/B)*b
//     for (int i = 0; i < 4; ++i) {
//         bStar[k][i] = (-u[k][k-1]) * bStar[k][i] + (B[k]/B_tmp)*b[i];
//     }

//     // B_k = B_{k-1}*B_k/B
//     B[k] = B[k-1]*B[k]/B_tmp;

//     // B_{k-1} = B
//     B[k-1] = B_tmp;

//     for (int i = k+1; i < kmax; ++i) {
//         double t = u[i][k];
//         // u_{i,k} = u_{i,k-1} - u*t
//         u[i][k] = u[i][k-1] - u_tmp*t;

//         // u_{i,k-1} = t + u_{k,k-1}*u_{i,k}
//         u[i][k-1] = t + u[k][k-1] * u[i][k];
//     }
// }

// static void dotproduct_row_dbl(int* mul, const int m1[4][4], const int m2[4][4], int m1j, int m2j) {
//     *mul = 0;
//     for (int i = 0; i < 4; ++i) {
//         *mul = *mul + m1[m1j][i]*m2[m2j][i];
//     }
// }

// static void dotproduct_zr_row_dbl(double* mul, const int m1[4][4], const double m2[4][4], int m1j, int m2j) {
//     *mul = 0;
//     for (int i = 0; i < 4; ++i) {
//         *mul = *mul + (double)m1[m1j][i]*m2[m2j][i];
//     }
// }

// static void dotproduct_rr_row_dbl(double* mul, const double m1[4][4], const double m2[4][4], int m1j, int m2j) {
//     *mul = 0;
//     for (int i = 0; i < 4; ++i) {
//         *mul = *mul + m1[m1j][i]*m2[m2j][i];
//     }
// }

// static void mul_row_dbl(double mul[4][4], const double* a, const double m[4][4], int j) {
//     for (int i = 0; i < 4; ++i) {
//         mul[j][i] = *a * m[j][i];
//     }
// }

// static void add_row_dbl(int add[4][4], const int a[4][4], const int b[4][4], int j, int aj, int bj) {
//     for (int i = 0; i < 4; ++i) {
//         add[j][i] = a[aj][i] +b[bj][i];
//     }
// }

// static void sub_row_dbl(double add[4][4], double a[4][4], double b[4][4], int j, int aj, int bj) {
//     for (int i = 0; i < 4; ++i) {
//         add[j][i] = a[aj][i] -b[bj][i];
//     }
// }

// /// @brief LLL reduction on 4-dimensional lattice
// /// Implements Algorithm 2.6.3 from Henri Cohen's "A Course in Computational Algebraic Number Theory"
// /// @param red 
// /// @param lattice 
// /// @return 
// int quat_lattice_lll_dbl(int red[4][4], const int lattice[4][4]) {
//     int ret = 0;
//     int basis[4][4];
//     double bStar[4][4];
//     double bStar_tmp[4][4];

//     double u[4][4];
//     int H[4][4]; // -> I_4
//     double B[4] = { 0 };

//     for (int i = 0; i < 4; ++i) {
//         for (int j = 0; j < 4; ++j) {
//             basis[i][j] = lattice[i][j];
//             u[i][j] = 0;
//             bStar[i][j] = 0;
//             bStar_tmp[i][j] = 0;
//             if (i == j)
//                 H[i][j] = 1;
//             else
//                 H[i][j] = 0;
//         }
//     }
//     int k = 1, kmax = 0;
//     // bStar_1 = b_1
//     for (int i = 0; i < 4; ++i)
//         bStar[0][i] = basis[0][i];

//     // B_1 = b_1 * b_1
//     int btmp = 0;
//     dotproduct_row_dbl(&btmp, basis, basis, 0, 0);
//     B[0] = (double) btmp;

//     while (k < 4) {
//         // Step 2: Incremental Gram-Schmidt
//         // if (k <= kmax) -> we can omit..
//         if (k > kmax) {
//             kmax = k;
//             for (int i = 0; i < 4; ++i) {

//                 bStar[k][i] = (double)basis[k][i];
//             }
//             for (int j = 0; j <= k-1; ++j) {
//                 // bStar_k = b_k -> already done initially -> todo: check if that's ok
//                 // nop

//                 // u_{k,j} = b_k*bSTAR_j/B_j
//                 dotproduct_zr_row_dbl(&u[k][j], basis, bStar, k, j);
//                 u[k][j] = u[k][j] / B[j];

//                 // bStar_k = bStar_k - u_{k,j}*bStar_j
//                 mul_row_dbl(bStar_tmp, &u[k][j], bStar, j);
//                 sub_row_dbl(bStar, bStar, bStar_tmp, k, k, j);
//             }
//             // B_k = bStar_k*bStar_k
//             dotproduct_rr_row_dbl(&B[k], bStar, bStar, k, k);
//             if (B[k] == 0) {
//                 // b_i did not form a basis, terminate with error
//                 ret = -1;
//                 goto err;
//             }
//         }

//         while(1) {
//             // Step 3: Test LLL condition
//             RED_dbl(basis, u, H, k, k - 1);
//             // If B_k < (0.75 - u_{k,k-1}^2)*B_{k-1}
//             double temp = (0.75 - u[k][k-1]*u[k][k-1])*B[k-1];
//             if (B[k] < temp) {
//                 SWAP_dbl(basis, u, H, B, bStar, k, kmax);
//                 k = (k - 1 > 1 ? k - 1 : 1);
//             } else {
//                 for (int l = k - 2; l >= 0; --l) { // TODO: or should it be k - 1 ? or k >= 2?
//                     RED_dbl(basis, u, H, k, l);
//                 }
//                 k++;
//                 break;
//             }
//         }
//     }

//     for (int i = 0; i < 4; ++i) {
//         for (int j = 0; j < 4; ++j) {
//             red[i][j] = basis[i][j];
//         }
//     }

// err:
//     return ret;
// }
// #endif