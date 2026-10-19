// Internal numerical utilities, GPL-2.0.
#pragma once
#include "bnsim/simulation.hpp"
#include <algorithm>
#include <exception>
#include <limits>
#include <stdexcept>
#include <thread>

namespace BNSim::detail {
// Each item writes only its own result; exceptions are propagated in item order.
template<class F> void parallel_for(std::size_t count, std::size_t workers, F f) {
    if (count == 0) return;
    workers = std::min(count, workers);
    if (workers == 1) { for (std::size_t i = 0; i < count; ++i) f(i); return; }
    std::vector<std::exception_ptr> errors(count);
    std::vector<std::thread> threads;
    threads.reserve(workers);
    try {
        for (std::size_t w = 0; w < workers; ++w)
            threads.emplace_back([&, w] {
                for (std::size_t i = count * w / workers; i < count * (w + 1) / workers; ++i)
                    try { f(i); } catch (...) { errors[i] = std::current_exception(); }
            });
    } catch (...) {
        for (auto& t : threads) t.join();
        throw;
    }
    for (auto& t : threads) t.join();
    for (const auto& e : errors) if (e) std::rethrow_exception(e);
}

struct OdeStats { std::size_t accepted = 0, rejected = 0, startup = 0; };
// Dormand-Prince 5(4), maximum scaled component error. Every stage and accepted
// endpoint must be finite and nonnegative. No concentration clipping. A
// step-doubled SSPRK3 startup handles zero-state reaction cascades for which
// non-SSP Dormand-Prince internal stages are negative at every step size.
template<class RHS> OdeStats integrate(std::vector<double>& y, double duration,
    const std::vector<double>& atol, const IntegrationOptions& options, RHS rhs) {
    const auto n = y.size();
    std::array<std::vector<double>, 7> k;
    for (auto& a : k) a.resize(n);
    std::vector<double> stage(n), fifth(n), fourth(n);
    auto valid = [](const auto& v) {
        return std::all_of(v.begin(), v.end(), [](double x) { return std::isfinite(x) && x >= 0; });
    };
    auto ssp = [&](const std::vector<double>& initial, double h, std::vector<double>& result) {
        std::vector<double> v=initial, derivative(n);
        for (int s=0; s<3; ++s) {
            std::fill(derivative.begin(),derivative.end(),0);
            rhs(v,derivative);
            for (std::size_t i=0;i<n;++i) {
                const double euler=v[i]+h*derivative[i];
                v[i]=s==0 ? euler : initial[i]+(s==1 ? .25 : 2./3)*(euler-initial[i]);
            }
            if (!valid(v)) return false;
        }
        result=std::move(v); return true;
    };
    if (!valid(y) || atol.size() != n) throw std::runtime_error("ODE invalid initial state");
    static constexpr double a[7][7] = {
        {}, {1./5}, {3./40, 9./40}, {44./45, -56./15, 32./9},
        {19372./6561, -25360./2187, 64448./6561, -212./729},
        {9017./3168, -355./33, 46732./5247, 49./176, -5103./18656},
        {35./384, 0, 500./1113, 125./192, -2187./6784, 11./84}};
    static constexpr double b4[7] = {5179./57600, 0, 7571./16695, 393./640, -92097./339200, 187./2100, 1./40};
    double time = 0, h = duration;
    OdeStats stats;
    for (std::size_t attempt = 0; time < duration && attempt < options.max_attempts; ++attempt) {
        h = std::min(h, duration - time);
        if (!(h > 0) || time + h == time) throw std::runtime_error("ODE timestep underflow");
        bool ok = true;
        std::size_t bad_stage=0, bad_component=0;
        double bad_value=0;
        for (std::size_t s = 0; s < 7; ++s) {
            for (std::size_t i = 0; i < n; ++i) {
                double sum = 0;
                for (std::size_t j = 0; j < s; ++j) sum += a[s][j] * k[j][i];
                stage[i] = y[i] + h * sum;
            }
            if (!valid(stage)) {
                ok = false; bad_stage=s;
                for (std::size_t i=0;i<n;++i) if (!std::isfinite(stage[i]) || stage[i]<0) { bad_component=i; bad_value=stage[i]; break; }
                break;
            }
            std::fill(k[s].begin(), k[s].end(), 0);
            rhs(stage, k[s]);
            if (!std::all_of(k[s].begin(), k[s].end(), [](double x) { return std::isfinite(x); })) { ok = false; break; }
        }
        double error = 0;
        bool startup=false;
        if (!ok && bad_value<0 && y[bad_component]==0) {
            // Algebraically negative internal stages in a reaction chain that
            // starts at zero cannot be cured by indefinite step halving.
            std::vector<double> full,half,twice;
            ++stats.rejected;
            if (ssp(y,h,full) && ssp(y,h/2,half) && ssp(half,h/2,twice)) {
                error=0;
                for (std::size_t i=0;i<n;++i)
                    error=std::max(error,std::abs(twice[i]-full[i])/(7*(atol[i]+options.relative_tolerance*std::max(y[i],twice[i]))));
                if (std::isfinite(error) && error<=1) { fifth=std::move(twice); ok=true; startup=true; }
            }
        }
        if (ok && !startup) {
            fifth = stage; // stage 7 is the fifth order endpoint
            for (std::size_t i = 0; i < n; ++i) {
                double sum = 0;
                for (std::size_t j = 0; j < 7; ++j) sum += b4[j] * k[j][i];
                fourth[i] = y[i] + h * sum;
                const double scale = atol[i] + options.relative_tolerance * std::max(y[i], fifth[i]);
                const double component_error=std::abs(fifth[i]-fourth[i])/scale;
                if (!std::isfinite(component_error)) { error=std::numeric_limits<double>::infinity(); break; }
                error = std::max(error,component_error);
            }
            ok = valid(fifth) && std::isfinite(error) && error <= 1;
        }
        if (ok) {
            y.swap(fifth);
            time = h == duration - time ? duration : time + h;
            ++stats.accepted;
            if (startup) ++stats.startup;
            h *= error == 0 ? 5 : std::clamp(.9 * std::pow(error, startup ? -.25 : -.2), .2, 5.);
        } else {
            ++stats.rejected;
            h *= error > 1 && std::isfinite(error) ? std::clamp(.9 * std::pow(error, -.2), .1, .5) : .5;
            if (h < options.minimum_step) {
                throw std::runtime_error("ODE cannot meet positivity/error tolerance above minimum step; stage="+std::to_string(bad_stage)+" component="+std::to_string(bad_component)+" value="+std::to_string(bad_value));
            }
        }
    }
    if (time < duration) throw std::runtime_error("ODE exceeded maximum attempts");
    return stats;
}
} // namespace BNSim::detail
