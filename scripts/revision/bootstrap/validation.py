"""A1 raw-scale validation method core; no automatic execution or file IO."""
import numpy as np
from aligned_dgm import generate
from aligned_geometry import contributions, aggregate, NAMES
from aligned_bootstrap import bootstrap_draws, infer
from oracle import oracle_band_classification


def validate_replication(config, cell_id, replication, population_endpoints, oracle_mcse,
                         B=None):
    """One prespecified replication against independently supplied population ratios.

    The frozen policies are synthetic pilot policies, never clinical models.
    Oracle arrays have shape (18, 2). K3 coverage contains all three contrasts.
    Calls perform computation only when explicitly invoked by a researcher.
    """
    cell = next(c for c in config['cells'] if c['id'] == cell_id)
    cfg = {**config['cell_defaults'], **cell}
    policy = config['policies'][cfg['policy_group']]
    plan = config['plans']['extended']
    if B is None:
        B = (plan['anchor_B'] if cell_id in plan['anchor_cell_ids'] and
             replication < plan['anchor_MC'] else plan['B'])
    data = generate(cfg, policy, np.random.SeedSequence(
        [config['seeds']['MC'], cfg['index'], replication]))
    num, den = contributions(data)
    point, sample_den = aggregate(num, den, np.ones(len(data)))
    point = point[0]
    draws, denominators = bootstrap_draws(
        num, den, data.hospitals, B,
        np.random.default_rng(np.random.SeedSequence(
            [config['seeds']['bootstrap'], cfg['index'], replication])),
        cfg['scheme'], config['inference']['batch_size'])
    output = []
    for report_B in sorted({min(plan['B'], B), B}):
        for mode in ('clinical_shared', 'target_specific'):
            rows = infer(point, draws[:report_B], mode,
                         config['inference']['minimum_valid_fraction'])
            outer = np.array([r['outer'] for r in rows])
            coverage = oracle_band_classification(outer, population_endpoints, oracle_mcse)
            for row in rows:
                j = row['index']
                output.append(dict(**row, target=NAMES[j], validity=mode,
                    sample_endpoints=point[j], sample_denominator=sample_den[0, j],
                    zero_denominator_draws=int((denominators[:report_B, j] == 0).sum()),
                    coverage=bool(coverage['nominal'][j]) if row['status'] == 'OK' else None,
                    coverage_definite=bool(coverage['definite'][j]) if row['status'] == 'OK' else None,
                    coverage_possible=bool(coverage['possible'][j]) if row['status'] == 'OK' else None))
            defined = all(r['status'] == 'OK' for r in rows[15:])
            output.append(dict(target='JOINT_A-B_A-C_B-C', validity=mode, B_requested=report_B,
                status='OK' if defined else 'UNDEFINED_OR_INVALID',
                coverage=bool(coverage['nominal'][15:].all()) if defined else None,
                coverage_definite=bool(coverage['definite'][15:].all()) if defined else None,
                coverage_possible=bool(coverage['possible'][15:].all()) if defined else None))
    return output
