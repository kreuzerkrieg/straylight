#pragma once

#include "Matrices.h"

void forward_model_single_psf_dual_resolution(Matrix2D& output_image, const Matrix2D& input_image, int channel, bool logStages = false);
