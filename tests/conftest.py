import os

# BEFORE numpy is imported (pytest loads this file before any test module): one BLAS thread.  The
# matrices of the suite are small and dense; with a multi-threaded BLAS the suite ran hundreds of
# times slower on a many-core machine.  A value set by the caller is kept.
for _var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_var, "1")

from pathlib import Path  # noqa: E402

import pytest  # noqa: E402

FIXTURES = Path(__file__).parent / "fixtures"


@pytest.fixture(autouse=True)
def _no_user_shapes_db(monkeypatch, tmp_path):
    """Tests never see a real user-supplied shapes database."""
    monkeypatch.setenv("BENCHFEA_SHAPES_DB", str(tmp_path / "absent-shapes-db.csv"))
    monkeypatch.delenv("BENCHFEA_SHAPES_COLUMNS", raising=False)
    monkeypatch.delenv("BENCHFEA_SHAPES_LENGTH_UNIT", raising=False)


@pytest.fixture
def fixtures():
    return FIXTURES


# ------------------------------------------------------------------------------------ slow tests
# ``pytest -m "not slow"`` is the contributor's run (about 2 minutes); plain ``pytest`` (CI, releases)
# runs everything.  The marker is given HERE, from measured durations, so that the list is in one place:
#
#   SLOW_FUNCTIONS   "file::function": every parametrisation of the test function
#   SLOW_IDS         "file::function[id]": one parametrisation
#   SLOW_GENERATED   ids of the generated models of tests/test_generators.py that cost seconds to make; every
#                    test parametrised with one of them is slow (the models are generated once per session
#                    and shared, so the first test that needs one pays for it)
#
# A test can also carry ``@pytest.mark.slow`` itself.  Nothing is skipped or weakened by the marker.  An
# entry that matches no collected test is reported as a warning (a renamed test must not drop out silently).
SLOW_FUNCTIONS = {
    "tests/test_dataset.py::test_cli_bench_make_dataset_and_run_agent",
    "tests/test_dataset.py::test_committed_tasks_are_solved_by_their_reference_answers",
    "tests/test_env.py::test_held_out_instances_are_refused_unless_the_environment_is_built_for_them",
    "tests/test_env.py::test_an_instance_without_a_verified_reference_is_refused_with_the_log",
    "tests/test_env.py::test_oracle_agent_scores_one_and_baseline_does_not",
    "tests/test_env.py::test_check_env",
    "tests/test_tasks.py::test_design_only_case_is_quarantined_not_served",
    "tests/test_tasks.py::test_design_task_and_reference_search_with_the_coverage_test_device",
    "tests/test_scoring.py::test_design_scoring_with_the_coverage_test_device",
    "tests/test_solve_roundoff_scales.py::test_buckling_result_says_when_the_lowest_mode_is_the_own_buckling_of_a_truss_member",
    "tests/test_solve_roundoff_scales.py::test_second_order_with_tension_only_bracing_under_gravity_is_solved_not_refused",
    "tests/test_solve_roundoff_scales.py::test_second_order_tension_only_result_against_opensees_and_pynite_with_the_same_active_set",
    "tests/test_solve_roundoff_scales.py::test_an_admissible_active_set_that_is_unstable_is_still_refused",
    "tests/test_crosscheck_estimate.py::test_estimated_error_is_not_below_the_true_error",
    "tests/test_crosscheck_estimate.py::test_rigid_diaphragm_in_second_order_pynite_is_unsupported_and_the_omitted_term_is_real",
    "tests/test_crosscheck_estimate.py::test_corrupted_engine_result_does_not_pass_tier_c",
    "tests/test_curveballs.py::test_targets_that_keep_the_model_cross_checkable_are_tried_first",
    "tests/test_curveballs.py::test_shared_sub_objects_of_a_python_document_are_not_edited_together",
    "tests/test_curveballs.py::test_every_generator_family_takes_curve_balls",
    "tests/test_generators.py::test_expected_cross_checkable_models_are_cross_checked",
    "tests/test_generators.py::test_opt_in_features_outside_the_cross_check_envelope",
    "tests/test_generators.py::test_units_and_vertical_axes_all_appear_and_describe_the_same_structure",
    "tests/test_generators.py::test_three_dimensional_families_are_measured_as_3d_and_plane_ones_as_planar",
    "tests/test_generators.py::test_family_specific_features",
    "tests/test_generators.py::test_tier_c_expectation_follows_the_features",
    "tests/test_verify_problems.py::test_library_of_derived_problems_and_its_quarantined_cases",
    # the four tests on a dataset built in the session (module fixture: four cells through the pipeline)
    "tests/test_dataset.py::test_build_writes_only_verified_tasks_and_logs_the_rest",
    "tests/test_dataset.py::test_build_is_deterministic_resumable_and_parallel",
    "tests/test_dataset.py::test_environment_can_serve_the_stored_dataset",
    "tests/test_dataset.py::test_cli_score",
    "tests/test_dataset.py::test_committed_task_is_reproduced_from_its_key",
    "tests/test_crosscheck_estimate.py::test_estimate_is_not_vacuous",
    "tests/test_crosscheck_estimate.py::test_solver_without_contraction_is_not_eligible",
    "tests/test_crosscheck_estimate.py::test_truss_member_with_an_axial_member_load_is_a_beam_column_in_the_foreign_models",
    "tests/test_splits.py::test_base_models_of_the_two_splits_are_different_documents",
    "tests/test_env.py::test_reset_without_options_draws_a_task_of_the_split",
    "tests/test_env.py::test_repair_in_the_environment_by_editing_the_working_copy",
    "tests/test_env.py::test_two_environments_with_the_same_seed_give_identical_transcripts",
    "tests/test_scoring.py::test_repair_is_resolved_through_the_pipeline_and_compared_with_the_intended_design",
    "tests/test_scoring.py::test_repair_compares_physical_values_whatever_the_units_of_the_submission",
    "tests/test_scoring.py::test_every_served_task_is_solvable_by_the_reference_answer",
    "tests/test_verify_design_coverage.py::test_registry_is_up_to_date_with_the_design_tests",
    "tests/test_generators.py::test_forced_features_appear_in_the_document_and_the_metadata",
    "tests/test_curveballs.py::test_apply_with_verification_accepts_only_a_verified_reference",
    "tests/test_curveballs.py::test_verification_route_is_recorded_and_degrading_targets_come_last",
    "tests/test_curveballs.py::test_seed_selects_the_target_and_nothing_else_decides",
    "tests/test_solve_stability.py::test_unsupported_features_are_refused_listing_every_one",
    "tests/test_verify_stations.py::test_axial_line_load_in_second_order_is_left_to_the_integration_and_passes",
    "tests/test_crosscheck_measured.py::test_opensees_timoshenko_element_mishandles_joint_offsets_and_ignores_pdelta",
    "tests/test_crosscheck_translators.py::test_second_order_portal_cross_check_between_solvers",
    "tests/independent/test_indep_second_order.py::test_second_order_with_varying_axial_force_reference_converges_to_the_engine",
    "tests/test_curveballs.py::test_every_curve_ball_is_applicable_to_a_generator_family",
    "tests/test_examples.py::test_example_05_only_the_stretched_diagonals_work",
    # users of the larger generated models of tests/test_curveballs.py (generated once, shared)
    "tests/test_curveballs.py::test_loaded_memberless_master_is_recorded_as_not_expected_to_verify",
    "tests/test_curveballs.py::test_second_order_off_only_where_the_amplification_is_material",
    "tests/test_curveballs.py::test_support_curve_balls_leave_the_plane_restraint_alone",
    "tests/test_curveballs.py::test_pinned_frame_is_a_mechanism_on_pinned_bases_and_a_cantilever_frame_on_fixed_ones",
    "tests/test_curveballs.py::test_constraint_made_flexible_replaces_the_constraint_by_members",
}
SLOW_IDS = {
    "tests/test_examples.py::test_example_runs_with_all_invariants_and_a_usable_design[05_braced_frame]",
    "tests/test_examples.py::test_documented_numbers_are_what_the_engine_gives[05_braced_frame]",
    "tests/independent/test_indep_one_sided.py::test_one_sided_members_against_the_tableau_reference"
    "[T02_space_tower_with_cross_bracing_y_up_kip_in-True]",
    "tests/test_env.py::test_reference_does_not_leak_before_the_episode_ends[change-repair]",
    "tests/test_tasks.py::test_description_states_every_value_of_the_model[braced_frame-1-2]",
    "tests/test_verify_stations.py::test_station_values_of_the_model_bank_are_reproduced[S07-source14]",
    "tests/test_verify_stations.py::test_station_values_of_the_model_bank_are_reproduced[S04-source11]",
    "tests/independent/test_indep_second_order.py::test_second_order_against_opensees_and_pynite"
    "[pynite-S04_skew_space_frame_y_up_kip_in-8]",
    "tests/independent/test_indep_second_order.py::test_second_order_against_opensees_and_pynite"
    "[opensees-S04_skew_space_frame_y_up_kip_in-32]",
} | {f"tests/test_curveballs.py::{test}[{name}]"
     for test in ("test_curve_ball_applies_record_describes_the_difference_and_effect_is_material",
                  "test_curve_ball_is_deterministic")
     for name in ("constraint_added", "constraint_removed", "constraint_made_flexible", "constraint_slaves_dropped",
                  "constraint_wrong_slave", "constraint_wrong_master", "second_order_off", "release_pinned_frame")}
SLOW_GENERATED = {
    "portal_frame-L2-s0", "portal_frame-L3-s1", "stair_frame-L2-s0", "stair_frame-L3-s1", "stair_frame-L4-s0",
    "braced_frame-L2-s2", "braced_frame-L3-s1", "multistorey_2d-L3-s0", "multistorey_3d-L3-s0",
    "multistorey_3d-L4-s0", "pipe_rack-L3-s0", "diaphragm_frame-L3-s2", "diaphragm_frame-L4-s2",
    "diaphragm_frame-L5-s0",
}


@pytest.hookimpl(tryfirst=True)
def pytest_collection_modifyitems(config, items):
    import warnings
    seen, files = set(), set()
    for item in items:
        nodeid = item.nodeid
        function, _, param = nodeid.partition("[")
        files.add(nodeid.split("::")[0])
        hit = None
        if function in SLOW_FUNCTIONS:
            hit = function
        elif nodeid in SLOW_IDS:
            hit = nodeid
        elif function.startswith("tests/test_generators.py::") and param.rstrip("]") in SLOW_GENERATED:
            hit = param.rstrip("]")
        if hit is not None:
            seen.add(hit)
            item.add_marker(pytest.mark.slow)
    if any("::" in str(a) for a in config.args):                 # single tests were asked for: no check
        return
    stale = sorted(e for e in SLOW_FUNCTIONS | SLOW_IDS if e.split("::")[0] in files and e not in seen)
    if "tests/test_generators.py" in files and not os.environ.get("BENCHFEA_GENERATOR_FULL"):
        stale += sorted(SLOW_GENERATED - seen)
    if stale:
        warnings.warn(pytest.PytestWarning("tests/conftest.py lists slow tests that were not collected: "
                                           + ", ".join(stale)))
