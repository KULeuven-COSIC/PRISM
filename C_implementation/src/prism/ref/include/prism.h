/** @file
 *
 * @authors Antonin Leroux
 *
 * @brief The protocols
 */

#ifndef PRISM_H
#define PRISM_H

#include <intbig.h>
#include <quaternion.h>
#include <klpt_constants.h>
#include <quaternion_data.h>
#include <torsion_constants.h>
#include <rng.h>
#include <endomorphism_action.h>
#include <encoded_sizes.h>
#include <fp_constants.h>

#include <stdio.h>
#include <klpt.h>
#include <id2iso.h>
#include <hd.h>
#include <dim2id2iso.h>


static const int exp_size = 219;

/** @defgroup primesign_primesign SQIsignHD protocols
 * @{
 */
/** @defgroup primesign_t Types for SQIsignHD protocols
 * @{
 */

/** @brief Type for the signature
 *
 * @typedef signature_t
 *
 * @struct signature
 *
 */

typedef struct signature
{
    ec_curve_t E_aux; /// the montgomery A coefficient for the auxilliary curve
    // ibz_mat_2x2_t mat_Baux_can_to_B_aux; /// the matrix of the desired basis
    ec_basis_t B_aux;
    int *hint_aux;
} signature_t;

/** @brief Type for the public keys
 *
 * @typedef public_key_t
 *
 * @struct public_key
 *
 */
typedef struct public_key
{
    ec_curve_t curve; /// the normalized A coefficient of the Montgomery curve
    int *hint_pk;
} public_key_t;

/** @brief Type for the secret keys
 *
 * @typedef secret_key_t
 *
 * @struct secret_key
 *
 */
typedef struct secret_key
{
    ec_basis_t canonical_basis; // the canonical basis of the public key curve
    ec_curve_t curve;           /// the public curve, but with little precomputations
    quat_left_ideal_t secret_ideal;
    ibz_mat_2x2_t mat_BAcan_to_BA0_two; /// mat_BA0_to_BAcan*BA0 = BAcan, where BAcan is the
                                        /// canonical basis of EA[2^e], and BA0 the image of the
                                        /// basis of E0[2^e] through the secret isogeny
} secret_key_t;

/** @}
 */

/*************************** Functions *****************************/

void protocols_keygen(public_key_t *pk, secret_key_t *sk);
int protocols_sign(signature_t *sig,
                   const public_key_t *pk,
                   secret_key_t *sk,
                   const unsigned char *m,
                   size_t l,
                   int verbose);
int protocols_verif(signature_t *sig, const public_key_t *pk, const unsigned char *m, size_t l);

void public_key_init(public_key_t *pk);
void public_key_finalize(public_key_t *pk);

void secret_key_init(secret_key_t *sk);
void secret_key_finalize(secret_key_t *sk);
void secret_sig_init(signature_t *sig);
void secret_sig_finalize(signature_t *sig);

void print_signature(const signature_t *sig);
void print_public_key(const public_key_t *pk);

/** @defgroup signature The signature protocol
 * @{
 */

/** @}
 */

/** @}
 */

#endif
