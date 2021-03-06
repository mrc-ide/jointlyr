// Generated by rstantools.  Do not edit by hand.

/*
    jointlyr is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    jointlyr is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with jointlyr.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef MODELS_HPP
#define MODELS_HPP
#define STAN__SERVICES__COMMAND_HPP
#include <rstan/rstaninc.hpp>
// Code generated by Stan version 2.21.0
#include <stan/model/model_header.hpp>
namespace model_rti0_bayesian_namespace {
using std::istream;
using std::string;
using std::stringstream;
using std::vector;
using stan::io::dump;
using stan::math::lgamma;
using stan::model::prob_grad;
using namespace stan::math;
static int current_statement_begin__;
stan::io::program_reader prog_reader__() {
    stan::io::program_reader reader;
    reader.add_event(0, 0, "start", "model_rti0_bayesian");
    reader.add_event(78, 76, "end", "model_rti0_bayesian");
    return reader;
}
template <typename T0__, typename T1__, typename T2__, typename T3__>
typename boost::math::tools::promote_args<T0__, T1__, T2__, T3__>::type
force_of_infection_obs(const Eigen::Matrix<T0__, Eigen::Dynamic, 1>& incid_init,
                           const Eigen::Matrix<T1__, Eigen::Dynamic, 1>& incid_obs,
                           const Eigen::Matrix<T2__, Eigen::Dynamic, 1>& omega,
                           const T3__& rt,
                           const int& back,
                           const int& si_trunc, std::ostream* pstream__) {
    typedef typename boost::math::tools::promote_args<T0__, T1__, T2__, T3__>::type local_scalar_t__;
    typedef local_scalar_t__ fun_return_scalar_t__;
    const static bool propto__ = true;
    (void) propto__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
    int current_statement_begin__ = -1;
    try {
        {
        current_statement_begin__ = 6;
        int T(0);
        (void) T;  // dummy to suppress unused var warning
        stan::math::fill(T, std::numeric_limits<int>::min());
        stan::math::assign(T,(rows(incid_init) + rows(incid_obs)));
        current_statement_begin__ = 7;
        validate_non_negative_index("incid_full", "T", T);
        Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> incid_full(T);
        stan::math::initialize(incid_full, DUMMY_VAR__);
        stan::math::fill(incid_full, DUMMY_VAR__);
        current_statement_begin__ = 8;
        local_scalar_t__ lambda(DUMMY_VAR__);
        (void) lambda;  // dummy to suppress unused var warning
        stan::math::initialize(lambda, DUMMY_VAR__);
        stan::math::fill(lambda, DUMMY_VAR__);
        current_statement_begin__ = 13;
        stan::model::assign(incid_full, 
                    stan::model::cons_list(stan::model::index_min_max(1, back), stan::model::nil_index_list()), 
                    incid_init, 
                    "assigning variable incid_full");
        current_statement_begin__ = 14;
        stan::model::assign(incid_full, 
                    stan::model::cons_list(stan::model::index_min_max((back + 1), rows(incid_full)), stan::model::nil_index_list()), 
                    incid_obs, 
                    "assigning variable incid_full");
        current_statement_begin__ = 15;
        stan::math::assign(lambda, dot_product(stan::model::rvalue(incid_full, stan::model::cons_list(stan::model::index_min_max((T - si_trunc), T), stan::model::nil_index_list()), "incid_full"), omega));
        current_statement_begin__ = 16;
        return stan::math::promote_scalar<fun_return_scalar_t__>(lambda);
        }
    } catch (const std::exception& e) {
        stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
        // Next line prevents compiler griping about no return
        throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
    }
}
struct force_of_infection_obs_functor__ {
    template <typename T0__, typename T1__, typename T2__, typename T3__>
        typename boost::math::tools::promote_args<T0__, T1__, T2__, T3__>::type
    operator()(const Eigen::Matrix<T0__, Eigen::Dynamic, 1>& incid_init,
                           const Eigen::Matrix<T1__, Eigen::Dynamic, 1>& incid_obs,
                           const Eigen::Matrix<T2__, Eigen::Dynamic, 1>& omega,
                           const T3__& rt,
                           const int& back,
                           const int& si_trunc, std::ostream* pstream__) const {
        return force_of_infection_obs(incid_init, incid_obs, omega, rt, back, si_trunc, pstream__);
    }
};
template <typename T0__, typename T1__, typename T2__>
Eigen::Matrix<typename boost::math::tools::promote_args<T0__, T1__, T2__>::type, Eigen::Dynamic, 1>
force_of_infection_unobs(const T0__& log_i0,
                             const Eigen::Matrix<T1__, Eigen::Dynamic, 1>& omega,
                             const T2__& rt,
                             const int& back,
                             const int& si_trunc, std::ostream* pstream__) {
    typedef typename boost::math::tools::promote_args<T0__, T1__, T2__>::type local_scalar_t__;
    typedef local_scalar_t__ fun_return_scalar_t__;
    const static bool propto__ = true;
    (void) propto__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
    int current_statement_begin__ = -1;
    try {
        {
        current_statement_begin__ = 22;
        validate_non_negative_index("incid_init", "back", back);
        Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> incid_init(back);
        stan::math::initialize(incid_init, DUMMY_VAR__);
        stan::math::fill(incid_init, DUMMY_VAR__);
        current_statement_begin__ = 23;
        int tmp(0);
        (void) tmp;  // dummy to suppress unused var warning
        stan::math::fill(tmp, std::numeric_limits<int>::min());
        current_statement_begin__ = 24;
        stan::math::assign(incid_init, rep_vector(0, back));
        current_statement_begin__ = 25;
        stan::model::assign(incid_init, 
                    stan::model::cons_list(stan::model::index_uni(1), stan::model::nil_index_list()), 
                    stan::math::exp(log_i0), 
                    "assigning variable incid_init");
        current_statement_begin__ = 26;
        for (int index = 2; index <= back; ++index) {
            current_statement_begin__ = 27;
            if (as_bool(logical_gt((index - si_trunc), 1))) {
                current_statement_begin__ = 27;
                stan::math::assign(tmp, (index - si_trunc));
            } else {
                current_statement_begin__ = 28;
                stan::math::assign(tmp, 1);
            }
            current_statement_begin__ = 29;
            stan::model::assign(incid_init, 
                        stan::model::cons_list(stan::model::index_uni(index), stan::model::nil_index_list()), 
                        (rt * dot_product(stan::model::rvalue(incid_init, stan::model::cons_list(stan::model::index_min_max(tmp, index), stan::model::nil_index_list()), "incid_init"), stan::model::rvalue(omega, stan::model::cons_list(stan::model::index_min_max(((si_trunc + 1) - (index - tmp)), (si_trunc + 1)), stan::model::nil_index_list()), "omega"))), 
                        "assigning variable incid_init");
        }
        current_statement_begin__ = 32;
        return stan::math::promote_scalar<fun_return_scalar_t__>(incid_init);
        }
    } catch (const std::exception& e) {
        stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
        // Next line prevents compiler griping about no return
        throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
    }
}
struct force_of_infection_unobs_functor__ {
    template <typename T0__, typename T1__, typename T2__>
        Eigen::Matrix<typename boost::math::tools::promote_args<T0__, T1__, T2__>::type, Eigen::Dynamic, 1>
    operator()(const T0__& log_i0,
                             const Eigen::Matrix<T1__, Eigen::Dynamic, 1>& omega,
                             const T2__& rt,
                             const int& back,
                             const int& si_trunc, std::ostream* pstream__) const {
        return force_of_infection_unobs(log_i0, omega, rt, back, si_trunc, pstream__);
    }
};
#include <stan_meta_header.hpp>
class model_rti0_bayesian
  : public stan::model::model_base_crtp<model_rti0_bayesian> {
private:
        int window;
        int window_back;
        vector_d incid;
        int si_trunc;
        vector_d omega;
        double log_incid_init;
public:
    model_rti0_bayesian(stan::io::var_context& context__,
        std::ostream* pstream__ = 0)
        : model_base_crtp(0) {
        ctor_body(context__, 0, pstream__);
    }
    model_rti0_bayesian(stan::io::var_context& context__,
        unsigned int random_seed__,
        std::ostream* pstream__ = 0)
        : model_base_crtp(0) {
        ctor_body(context__, random_seed__, pstream__);
    }
    void ctor_body(stan::io::var_context& context__,
                   unsigned int random_seed__,
                   std::ostream* pstream__) {
        typedef double local_scalar_t__;
        boost::ecuyer1988 base_rng__ =
          stan::services::util::create_rng(random_seed__, 0);
        (void) base_rng__;  // suppress unused var warning
        current_statement_begin__ = -1;
        static const char* function__ = "model_rti0_bayesian_namespace::model_rti0_bayesian";
        (void) function__;  // dummy to suppress unused var warning
        size_t pos__;
        (void) pos__;  // dummy to suppress unused var warning
        std::vector<int> vals_i__;
        std::vector<double> vals_r__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
        try {
            // initialize data block variables from context__
            current_statement_begin__ = 36;
            context__.validate_dims("data initialization", "window", "int", context__.to_vec());
            window = int(0);
            vals_i__ = context__.vals_i("window");
            pos__ = 0;
            window = vals_i__[pos__++];
            current_statement_begin__ = 38;
            context__.validate_dims("data initialization", "window_back", "int", context__.to_vec());
            window_back = int(0);
            vals_i__ = context__.vals_i("window_back");
            pos__ = 0;
            window_back = vals_i__[pos__++];
            current_statement_begin__ = 42;
            validate_non_negative_index("incid", "window", window);
            context__.validate_dims("data initialization", "incid", "vector_d", context__.to_vec(window));
            incid = Eigen::Matrix<double, Eigen::Dynamic, 1>(window);
            vals_r__ = context__.vals_r("incid");
            pos__ = 0;
            size_t incid_j_1_max__ = window;
            for (size_t j_1__ = 0; j_1__ < incid_j_1_max__; ++j_1__) {
                incid(j_1__) = vals_r__[pos__++];
            }
            current_statement_begin__ = 43;
            context__.validate_dims("data initialization", "si_trunc", "int", context__.to_vec());
            si_trunc = int(0);
            vals_i__ = context__.vals_i("si_trunc");
            pos__ = 0;
            si_trunc = vals_i__[pos__++];
            current_statement_begin__ = 46;
            validate_non_negative_index("omega", "(si_trunc + 1)", (si_trunc + 1));
            context__.validate_dims("data initialization", "omega", "vector_d", context__.to_vec((si_trunc + 1)));
            omega = Eigen::Matrix<double, Eigen::Dynamic, 1>((si_trunc + 1));
            vals_r__ = context__.vals_r("omega");
            pos__ = 0;
            size_t omega_j_1_max__ = (si_trunc + 1);
            for (size_t j_1__ = 0; j_1__ < omega_j_1_max__; ++j_1__) {
                omega(j_1__) = vals_r__[pos__++];
            }
            current_statement_begin__ = 48;
            context__.validate_dims("data initialization", "log_incid_init", "double", context__.to_vec());
            log_incid_init = double(0);
            vals_r__ = context__.vals_r("log_incid_init");
            pos__ = 0;
            log_incid_init = vals_r__[pos__++];
            // initialize transformed data variables
            // execute transformed data statements
            // validate transformed data
            // validate, set parameter ranges
            num_params_r__ = 0U;
            param_ranges_i__.clear();
            current_statement_begin__ = 52;
            num_params_r__ += 1;
            current_statement_begin__ = 53;
            num_params_r__ += 1;
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }
    ~model_rti0_bayesian() { }
    void transform_inits(const stan::io::var_context& context__,
                         std::vector<int>& params_i__,
                         std::vector<double>& params_r__,
                         std::ostream* pstream__) const {
        typedef double local_scalar_t__;
        stan::io::writer<double> writer__(params_r__, params_i__);
        size_t pos__;
        (void) pos__; // dummy call to supress warning
        std::vector<double> vals_r__;
        std::vector<int> vals_i__;
        current_statement_begin__ = 52;
        if (!(context__.contains_r("log_i0")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable log_i0 missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("log_i0");
        pos__ = 0U;
        context__.validate_dims("parameter initialization", "log_i0", "double", context__.to_vec());
        double log_i0(0);
        log_i0 = vals_r__[pos__++];
        try {
            writer__.scalar_lub_unconstrain((log_incid_init - ((stan::math::log(10) * 100) / 6)), (log_incid_init - ((stan::math::log(0.5) * 100) / 6)), log_i0);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable log_i0: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        current_statement_begin__ = 53;
        if (!(context__.contains_r("rt_est")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable rt_est missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("rt_est");
        pos__ = 0U;
        context__.validate_dims("parameter initialization", "rt_est", "double", context__.to_vec());
        double rt_est(0);
        rt_est = vals_r__[pos__++];
        try {
            writer__.scalar_lub_unconstrain(0.1, 10, rt_est);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable rt_est: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        params_r__ = writer__.data_r();
        params_i__ = writer__.data_i();
    }
    void transform_inits(const stan::io::var_context& context,
                         Eigen::Matrix<double, Eigen::Dynamic, 1>& params_r,
                         std::ostream* pstream__) const {
      std::vector<double> params_r_vec;
      std::vector<int> params_i_vec;
      transform_inits(context, params_i_vec, params_r_vec, pstream__);
      params_r.resize(params_r_vec.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r(i) = params_r_vec[i];
    }
    template <bool propto__, bool jacobian__, typename T__>
    T__ log_prob(std::vector<T__>& params_r__,
                 std::vector<int>& params_i__,
                 std::ostream* pstream__ = 0) const {
        typedef T__ local_scalar_t__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // dummy to suppress unused var warning
        T__ lp__(0.0);
        stan::math::accumulator<T__> lp_accum__;
        try {
            stan::io::reader<local_scalar_t__> in__(params_r__, params_i__);
            // model parameters
            current_statement_begin__ = 52;
            local_scalar_t__ log_i0;
            (void) log_i0;  // dummy to suppress unused var warning
            if (jacobian__)
                log_i0 = in__.scalar_lub_constrain((log_incid_init - ((stan::math::log(10) * 100) / 6)), (log_incid_init - ((stan::math::log(0.5) * 100) / 6)), lp__);
            else
                log_i0 = in__.scalar_lub_constrain((log_incid_init - ((stan::math::log(10) * 100) / 6)), (log_incid_init - ((stan::math::log(0.5) * 100) / 6)));
            current_statement_begin__ = 53;
            local_scalar_t__ rt_est;
            (void) rt_est;  // dummy to suppress unused var warning
            if (jacobian__)
                rt_est = in__.scalar_lub_constrain(0.1, 10, lp__);
            else
                rt_est = in__.scalar_lub_constrain(0.1, 10);
            // model body
            {
            current_statement_begin__ = 56;
            validate_non_negative_index("incid_est", "window_back", window_back);
            Eigen::Matrix<local_scalar_t__, Eigen::Dynamic, 1> incid_est(window_back);
            stan::math::initialize(incid_est, DUMMY_VAR__);
            stan::math::fill(incid_est, DUMMY_VAR__);
            current_statement_begin__ = 57;
            local_scalar_t__ lambda(DUMMY_VAR__);
            (void) lambda;  // dummy to suppress unused var warning
            stan::math::initialize(lambda, DUMMY_VAR__);
            stan::math::fill(lambda, DUMMY_VAR__);
            current_statement_begin__ = 58;
            stan::math::assign(incid_est, force_of_infection_unobs(log_i0, omega, rt_est, window_back, si_trunc, pstream__));
            current_statement_begin__ = 61;
            for (int index = 1; index <= window; ++index) {
                current_statement_begin__ = 62;
                stan::math::assign(lambda, force_of_infection_obs(incid_est, stan::model::rvalue(incid, stan::model::cons_list(stan::model::index_min_max(1, index), stan::model::nil_index_list()), "incid"), omega, rt_est, window_back, si_trunc, pstream__));
                current_statement_begin__ = 64;
                lp_accum__.add(((-(rt_est) * lambda) + (get_base1(incid, index, "incid", 1) * stan::math::log((rt_est * lambda)))));
            }
            }
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
        lp_accum__.add(lp__);
        return lp_accum__.sum();
    } // log_prob()
    template <bool propto, bool jacobian, typename T_>
    T_ log_prob(Eigen::Matrix<T_,Eigen::Dynamic,1>& params_r,
               std::ostream* pstream = 0) const {
      std::vector<T_> vec_params_r;
      vec_params_r.reserve(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        vec_params_r.push_back(params_r(i));
      std::vector<int> vec_params_i;
      return log_prob<propto,jacobian,T_>(vec_params_r, vec_params_i, pstream);
    }
    void get_param_names(std::vector<std::string>& names__) const {
        names__.resize(0);
        names__.push_back("log_i0");
        names__.push_back("rt_est");
        names__.push_back("incid_est");
    }
    void get_dims(std::vector<std::vector<size_t> >& dimss__) const {
        dimss__.resize(0);
        std::vector<size_t> dims__;
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dims__.push_back(((window_back + window) + 7));
        dimss__.push_back(dims__);
    }
    template <typename RNG>
    void write_array(RNG& base_rng__,
                     std::vector<double>& params_r__,
                     std::vector<int>& params_i__,
                     std::vector<double>& vars__,
                     bool include_tparams__ = true,
                     bool include_gqs__ = true,
                     std::ostream* pstream__ = 0) const {
        typedef double local_scalar_t__;
        vars__.resize(0);
        stan::io::reader<local_scalar_t__> in__(params_r__, params_i__);
        static const char* function__ = "model_rti0_bayesian_namespace::write_array";
        (void) function__;  // dummy to suppress unused var warning
        // read-transform, write parameters
        double log_i0 = in__.scalar_lub_constrain((log_incid_init - ((stan::math::log(10) * 100) / 6)), (log_incid_init - ((stan::math::log(0.5) * 100) / 6)));
        vars__.push_back(log_i0);
        double rt_est = in__.scalar_lub_constrain(0.1, 10);
        vars__.push_back(rt_est);
        double lp__ = 0.0;
        (void) lp__;  // dummy to suppress unused var warning
        stan::math::accumulator<double> lp_accum__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
        if (!include_tparams__ && !include_gqs__) return;
        try {
            if (!include_gqs__ && !include_tparams__) return;
            if (!include_gqs__) return;
            // declare and define generated quantities
            current_statement_begin__ = 70;
            validate_non_negative_index("incid_est", "((window_back + window) + 7)", ((window_back + window) + 7));
            Eigen::Matrix<double, Eigen::Dynamic, 1> incid_est(((window_back + window) + 7));
            stan::math::initialize(incid_est, DUMMY_VAR__);
            stan::math::fill(incid_est, DUMMY_VAR__);
            // generated quantities statements
            current_statement_begin__ = 72;
            stan::math::assign(incid_est, force_of_infection_unobs(log_i0, omega, rt_est, ((window_back + window) + 7), si_trunc, pstream__));
            // validate, write generated quantities
            current_statement_begin__ = 70;
            size_t incid_est_j_1_max__ = ((window_back + window) + 7);
            for (size_t j_1__ = 0; j_1__ < incid_est_j_1_max__; ++j_1__) {
                vars__.push_back(incid_est(j_1__));
            }
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }
    template <typename RNG>
    void write_array(RNG& base_rng,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& params_r,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& vars,
                     bool include_tparams = true,
                     bool include_gqs = true,
                     std::ostream* pstream = 0) const {
      std::vector<double> params_r_vec(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r_vec[i] = params_r(i);
      std::vector<double> vars_vec;
      std::vector<int> params_i_vec;
      write_array(base_rng, params_r_vec, params_i_vec, vars_vec, include_tparams, include_gqs, pstream);
      vars.resize(vars_vec.size());
      for (int i = 0; i < vars.size(); ++i)
        vars(i) = vars_vec[i];
    }
    std::string model_name() const {
        return "model_rti0_bayesian";
    }
    void constrained_param_names(std::vector<std::string>& param_names__,
                                 bool include_tparams__ = true,
                                 bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        param_name_stream__.str(std::string());
        param_name_stream__ << "log_i0";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "rt_est";
        param_names__.push_back(param_name_stream__.str());
        if (!include_gqs__ && !include_tparams__) return;
        if (include_tparams__) {
        }
        if (!include_gqs__) return;
        size_t incid_est_j_1_max__ = ((window_back + window) + 7);
        for (size_t j_1__ = 0; j_1__ < incid_est_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "incid_est" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
    }
    void unconstrained_param_names(std::vector<std::string>& param_names__,
                                   bool include_tparams__ = true,
                                   bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        param_name_stream__.str(std::string());
        param_name_stream__ << "log_i0";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "rt_est";
        param_names__.push_back(param_name_stream__.str());
        if (!include_gqs__ && !include_tparams__) return;
        if (include_tparams__) {
        }
        if (!include_gqs__) return;
        size_t incid_est_j_1_max__ = ((window_back + window) + 7);
        for (size_t j_1__ = 0; j_1__ < incid_est_j_1_max__; ++j_1__) {
            param_name_stream__.str(std::string());
            param_name_stream__ << "incid_est" << '.' << j_1__ + 1;
            param_names__.push_back(param_name_stream__.str());
        }
    }
}; // model
}  // namespace
typedef model_rti0_bayesian_namespace::model_rti0_bayesian stan_model;
#ifndef USING_R
stan::model::model_base& new_model(
        stan::io::var_context& data_context,
        unsigned int seed,
        std::ostream* msg_stream) {
  stan_model* m = new stan_model(data_context, seed, msg_stream);
  return *m;
}
#endif
#endif
