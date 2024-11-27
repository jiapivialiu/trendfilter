# admm control lists check args

    Code
      admm_control_list(max_iter = 0L)
    Condition
      Error in `admm_control_list()`:
      ! Assertion on 'max_iter' failed: Element 1 is not >= 1.

---

    Code
      admm_control_list(max_iter = 1.1)
    Condition
      Error in `admm_control_list()`:
      ! Assertion on 'max_iter' failed: Must be of type 'integerish', but element 1 is not close to an integer.

---

    Code
      admm_control_list(max_iter = 100:200)
    Condition
      Error in `admm_control_list()`:
      ! Assertion on 'max_iter' failed: Must have length 1, but has length 101.

---

    Code
      admm_control_list(rho_scale = 0)
    Condition
      Error in `admm_control_list()`:
      ! Assertion on 'rho_scale' failed: Element 1 is not >= 2.22045e-16.

---

    Code
      admm_control_list(rho_scale = Inf)
    Condition
      Error in `admm_control_list()`:
      ! Assertion on 'rho_scale' failed: Must be finite.

---

    Code
      admm_control_list(rho_scale = c(1e-06, 1e-05))
    Condition
      Error in `admm_control_list()`:
      ! Assertion on 'rho_scale' failed: Must have length 1, but has length 2.

---

    Code
      admm_control_list(tolerance = 0)
    Condition
      Error in `admm_control_list()`:
      ! Assertion on 'tolerance' failed: Element 1 is not >= 2.22045e-16.

---

    Code
      admm_control_list(tolerance = Inf)
    Condition
      Error in `admm_control_list()`:
      ! Assertion on 'tolerance' failed: Must be finite.

---

    Code
      admm_control_list(tolerance = c(1e-06, 1e-05))
    Condition
      Error in `admm_control_list()`:
      ! Assertion on 'tolerance' failed: Must have length 1, but has length 2.

---

    Code
      admm_control_list(unknown = 1)
    Condition
      Error in `admm_control_list()`:
      ! `...` must be empty.
      x Problematic argument:
      * unknown = 1

# trendfilter control lists check args

    Code
      trendfilter_control_list(obj_tol = 0)
    Condition
      Error in `trendfilter_control_list()`:
      ! Assertion on 'obj_tol' failed: Element 1 is not >= 2.22045e-16.

---

    Code
      trendfilter_control_list(obj_tol = Inf)
    Condition
      Error in `trendfilter_control_list()`:
      ! Assertion on 'obj_tol' failed: Must be finite.

---

    Code
      trendfilter_control_list(obj_tol = c(1e-06, 1e-05))
    Condition
      Error in `trendfilter_control_list()`:
      ! Assertion on 'obj_tol' failed: Must have length 1, but has length 2.

---

    Code
      trendfilter_control_list(x_cond = 0)
    Condition
      Error in `trendfilter_control_list()`:
      ! Assertion on 'x_cond' failed: Element 1 is not >= 1.

---

    Code
      trendfilter_control_list(x_cond = Inf)
    Condition
      Error in `trendfilter_control_list()`:
      ! Assertion on 'x_cond' failed: Must be finite.

---

    Code
      trendfilter_control_list(x_cond = c(5, 10))
    Condition
      Error in `trendfilter_control_list()`:
      ! Assertion on 'x_cond' failed: Must have length 1, but has length 2.

---

    Code
      trendfilter_control_list(admm_control = list())
    Condition
      Error in `trendfilter_control_list()`:
      ! Assertion on 'admm_control' failed: Must inherit from class 'admm_control', but has class 'list'.

---

    Code
      trendfilter_control_list(unknown = 1)
    Condition
      Error in `trendfilter_control_list()`:
      ! `...` must be empty.
      x Problematic argument:
      * unknown = 1

