"""Given one atom, compute the Laplace approximation to the posterior at 
that atom, where the observed data is defined to be the noise-free data generated 
by the chosen atom.

One can show that, at the chosen atom, the gradient of the log likelihood is zero and 
the exact Hessian of the log likelihood is

    d^2 log p(data | theta) / d theta^2 = -J.T @ J / sigma**2,

where ``J`` is the analytic Jacobian of the forward map with respect to
``theta = (x0, depth, R, amplitude)``. The exponential radius prior gives the 
log prior a nonzero gradient and zero Hessian. The Laplace approximation to 
the posterior therefore has covariance equal to the inverse of the Hessian of 
the log likelihood and has mean 

    mean = theta + covariance @ gradient_log_posterior = theta + covariance @ gradient_log_prior.

The atom, noise level, scene path, and output path are all chosen directly in ``main``.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from bayesian_gpr import Atom, Scene
from bayesian_gpr.forward import render_atom


BASIN_DIR = Path(__file__).resolve().parent
DATA_DIR = BASIN_DIR / "data"

PARAMETER_NAMES = ("x0", "depth", "R", "amplitude")


@dataclass(frozen=True)
class LaplaceApproximation:
    """Store a one-atom Laplace approximation and supporting quantities.

    Attributes:
        atom: Chosen representative atom.
        jacobian: Analytic forward-map Jacobian with shape ``(Nt * Nx, 4)``.
        gradient_log_posterior: Gradient of the log posterior at the representative atom.
        mean: Mean of the Gaussian Laplace approximation.
        precision: Negative Hessian of the log posterior.
        covariance: Moore--Penrose pseudoinverse of the precision matrix.
        precision_eigenvalues: Eigenvalues of the precision matrix.
        condition_number: Condition number of the precision matrix.
    """

    atom: Atom
    jacobian: np.ndarray
    gradient_log_posterior: np.ndarray
    mean: np.ndarray
    precision: np.ndarray
    covariance: np.ndarray
    precision_eigenvalues: np.ndarray
    condition_number: float


def ricker_and_derivative(
    tau: np.ndarray,
    f0: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Return the Ricker wavelet and its derivative with respect to time offset.

    Args:
        tau: Time offsets in seconds.
        f0: Centre frequency in Hz.

    Returns:
        Ricker-wavelet values and their derivatives with respect to ``tau``.
    """
    c_squared = (np.pi * f0) ** 2
    u = c_squared * tau**2
    exp_minus_u = np.exp(-u)

    wavelet = (1.0 - 2.0 * u) * exp_minus_u
    wavelet_derivative = (
        2.0 * c_squared * tau * (2.0 * u - 3.0) * exp_minus_u
    )
    return wavelet, wavelet_derivative


def analytic_atom_jacobian(
    atom: Atom,
    scene: Scene,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Analytically compute the Jacobian of the forward map with 
    respect to the atom parameters, evaluated at the given atom.

    The parameter order is ``(x0, depth, R, amplitude)``. Each Jacobian 
    column is the vectorised derivative of the forward map output with 
    respect to one atom parameter.

    Args:
        atom: Atom at which to evaluate the Jacobian of the forward map.
        scene: Scene supplying the grid, wavelet, and soil parameters.

    Returns:
        Jacobian of the forward map with shape ``(Nt * Nx, 4)``.
    """
    x_axis = np.asarray(scene.grid.x_axis, dtype=float)
    t_axis = np.asarray(scene.grid.t_axis, dtype=float)
    velocity = float(scene.soil.v)
    alpha = float(scene.soil.alpha)
    f0 = float(scene.wavelet.f0)

    effective_depth = float(atom.depth + atom.R)
    if effective_depth <= 0.0:
        raise ValueError("Atom depth + radius must be positive")
    if velocity <= 0.0:
        raise ValueError("Soil velocity must be positive")

    horizontal_offset = x_axis - float(atom.x0)
    slant_range = np.sqrt(effective_depth**2 + horizontal_offset**2)

    spread = np.sqrt(effective_depth / slant_range)
    absorption = np.exp(-2.0 * alpha * slant_range)
    attenuation = spread * absorption

    # This is the same travel-time locus used by render_atom.
    travel_time = (
        (2.0 / velocity) * (slant_range - effective_depth)
        + (2.0 * atom.depth / velocity)
    )
    tau = t_axis[:, None] - travel_time[None, :]

    wavelet, wavelet_derivative = ricker_and_derivative(tau, f0)
    amplitude = float(atom.amplitude)

    # Derivatives of log(B), where
    # B = sqrt(effective_depth / slant_range) * exp(-2 alpha slant_range).
    dlog_attenuation_dx0 = horizontal_offset * (
        0.5 / slant_range**2 + 2.0 * alpha / slant_range
    )
    dlog_attenuation_deffective_depth = (
        0.5 / effective_depth
        - 0.5 * effective_depth / slant_range**2
        - 2.0 * alpha * effective_depth / slant_range
    )

    # Since tau = t - travel_time, these are the derivatives of tau.
    dtau_dx0 = 2.0 * horizontal_offset / (velocity * slant_range)
    dtau_ddepth = -2.0 * effective_depth / (velocity * slant_range)
    dtau_dR = (2.0 / velocity) * (
        1.0 - effective_depth / slant_range
    )

    common = amplitude * attenuation[None, :]

    dA_dx0 = common * (
        wavelet * dlog_attenuation_dx0[None, :]
        + wavelet_derivative * dtau_dx0[None, :]
    )
    dA_ddepth = common * (
        wavelet * dlog_attenuation_deffective_depth[None, :]
        + wavelet_derivative * dtau_ddepth[None, :]
    )
    dA_dR = common * (
        wavelet * dlog_attenuation_deffective_depth[None, :]
        + wavelet_derivative * dtau_dR[None, :]
    )
    dA_damplitude = attenuation[None, :] * wavelet

    derivative_images = np.stack(
        (dA_dx0, dA_ddepth, dA_dR, dA_damplitude),
        axis=0,
    )
    jacobian = derivative_images.reshape(4, -1).T

    return jacobian


def compute_laplace_approximation(
    sigma: float,
    atom: Atom,
    scene: Scene,
    *,
    lambda_r: float = 10.0,
) -> LaplaceApproximation:
    """Compute the one-atom Jacobian, Hessian, and covariance. TODO: CHANGE DOCSTRING

    The data are always generated from the same atom at which the derivatives
    are evaluated:

        data = A(atom).

    Consequently, the residual is zero and the exact Hessian of the Gaussian
    log likelihood is ``-J.T @ J / sigma**2``. The uniform priors add no 
    gradient or Hessian curvature. The exponential radius prior has log
    density ``constant - lambda_r * R``. It therefore contributes ``-lambda_r``
    to the radius component of the gradient, but contributes no Hessian curvature.

    Completing the square gives the recentered local Gaussian

        N(mean, covariance),

    where ``mean = theta + covariance @ gradient_log_posterior``.

    Args:
        sigma: Standard deviation of the additive Gaussian noise model.
        atom: Representative atom defining both the synthetic data and the
            point at which curvature is calculated.
        scene: Scene supplying the grid, wavelet, and soil parameters.
        lambda_r: Rate of the exponential prior on the radius.

    Returns:
        one-atom Laplace approximation and supporting quantities.
    """
    if sigma <= 0.0 or not np.isfinite(sigma):
        raise ValueError("sigma must be finite and positive")
    if lambda_r < 0.0 or not np.isfinite(lambda_r):
        raise ValueError("lambda_r must be finite and nonnegative")

    gradient_log_prior = np.zeros(len(PARAMETER_NAMES), dtype=float)
    gradient_log_prior[2] = -lambda_r
    gradient_log_posterior = gradient_log_prior

    jacobian = analytic_atom_jacobian(atom, scene)
    hessian_log_likelihood = -(jacobian.T @ jacobian) / sigma**2
    precision = -hessian_log_likelihood
    eigenvalues = np.linalg.eigvalsh(precision)
    condition_number = float(np.linalg.cond(precision))

    # Use the pseudoinverse because a parameter direction may be locally
    # unidentifiable, in which case the precision matrix is singular.
    covariance = np.linalg.pinv(precision, hermitian=True)

    theta = np.array(
        [atom.x0, atom.depth, atom.R, atom.amplitude],
        dtype=float,
    )
    mean = theta + covariance @ gradient_log_posterior

    return LaplaceApproximation(
        atom = atom,
        jacobian = jacobian,
        gradient_log_posterior = gradient_log_posterior,
        mean = mean,
        precision = precision,
        covariance = covariance,
        precision_eigenvalues = eigenvalues,
        condition_number = condition_number,
    )


def main() -> None:
    """Choose an atom and compute its local Laplace approximation.

    Args:
        None.

    Returns:
        None.
    """
    # -------------------------------------------------------------------------
    # User choices
    # -------------------------------------------------------------------------
    scene_path = DATA_DIR / "scene.npz"
    output_path = BASIN_DIR / "atom_curvature.npz"

    sigma = 0.02
    lambda_r = 10.0

    atom = Atom(
        x0=25.0,
        depth=1.5,
        R=0.10,
        amplitude=0.8,
        eps_r=1.0,
    )

    # -------------------------------------------------------------------------
    # Laplace approximation
    # -------------------------------------------------------------------------
    scene = Scene.load(scene_path)
    result = compute_laplace_approximation(
        sigma,
        atom,
        scene,
        lambda_r=lambda_r,
    )

    print("Parameter order:", PARAMETER_NAMES)
    print("Atom:")
    print(np.array([atom.x0, atom.depth, atom.R, atom.amplitude]))
    print("\nGradient of log posterior:")
    print(result.gradient_log_posterior)
    print("\nMean:")
    print(result.mean)
    print("\nPrecision:")
    print(result.precision)
    print("\nCovariance (Moore--Penrose inverse of precision):")
    print(result.covariance)
    print("\nPrecision eigenvalues:")
    print(result.precision_eigenvalues)
    print(f"Condition number: {result.condition_number:.6g}")

    # -------------------------------------------------------------------------
    # Save results
    # -------------------------------------------------------------------------
    output_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        output_path,
        parameter_names=np.array(PARAMETER_NAMES),
        sigma=sigma,
        lambda_r=lambda_r,
        atom=np.array([atom.x0, atom.depth, atom.R, atom.amplitude]),
        jacobian=result.jacobian,
        gradient_log_posterior=result.gradient_log_posterior,
        mean=result.mean,
        precision=result.precision,
        covariance=result.covariance,
        precision_eigenvalues=result.precision_eigenvalues,
        condition_number=result.condition_number,
    )
    print(f"\nSaved results to {output_path.resolve()}")


if __name__ == "__main__":
    main()