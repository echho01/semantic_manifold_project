"""Finite-alphabet conditional-source experiment; beta uses bit distortion and bit information."""

import torch
from reference_graph import coordinates, fit_graph, attached
from link import (
    arguments,
    device,
    train_source,
    test_source,
    kl_table,
    information,
    report,
    save,
    LN2,
    post,
)


def solve(cost, beta, iters, tol):
    cost = cost.double()
    logk = -LN2 * beta * cost
    logm = torch.full_like(cost[0], -torch.log(cost.new_tensor(cost.shape[1])))
    for step in range(1, iters + 1):
        logz = torch.logsumexp(logm + logk, 1)
        lw = logm + logk - logz[:, None]
        W = lw.exp()
        D = (W * cost).sum(1).mean()
        upper = information(W) + beta * D
        loga = torch.logsumexp(logk - logz[:, None], 0) - torch.log(
            cost.new_tensor(len(cost))
        )
        lower = (-logz.mean() - loga.max()) / LN2
        gap = float((upper - lower).clamp_min(0))
        if gap <= tol or step == iters:
            break
        logm = (logm + loga).log_softmax(0)
    return logm.cpu(), dict(
        iterations=step,
        rate_bits=float(information(W)),
        distortion=float(D),
        upper=float(upper),
        lower=float(lower),
        gap=gap,
        tolerance_met=gap <= tol,
    )


def main():
    ap, tr, _ = arguments(__doc__)
    tr.add_argument(
        "--criterion",
        choices=("posterior_kl", "task_chord", "graph_geodesic"),
        default="posterior_kl",
    )
    tr.add_argument(
        "--eta",
        type=float,
        default=0.1,
        help="fixed reference weight in product task coordinates",
    )
    tr.add_argument(
        "--neighbors",
        type=int,
        default=12,
        help="training union kNN and out-of-sample attachment",
    )
    tr.add_argument("--beta", type=float, default=3.0)
    tr.add_argument(
        "--alphabet",
        choices=("objects", "augmented"),
        default="objects",
        help="augmented adds one free simplex output per anchor: midpoint with next anchor",
    )
    tr.add_argument("--iterations", type=int, default=500)
    tr.add_argument("--tolerance", type=float, default=1e-5)
    a = ap.parse_args()
    dev = device(a.device)
    if a.command == "train":
        if a.beta < 0 or a.iterations < 1 or a.tolerance <= 0:
            raise ValueError("invalid solver settings")
        ck, x, q = train_source(a, dev)
        b = post(ck["anchors"].to(dev), ck)
        if a.alphabet == "augmented":
            b = torch.cat((b, (b + b.roll(-1, 0)) / 2))
        q, b = q.double(), b.double()
        if a.eta <= 0 or (a.criterion != "posterior_kl" and a.alphabet != "objects"):
            raise ValueError(
                "eta must be positive; geometric costs need object-valued outputs"
            )
        ck.update(criterion=a.criterion, eta=a.eta)
        if a.criterion == "posterior_kl":
            cost = kl_table(q, b)
            unit = "posterior_KL_bits"
        else:
            u = coordinates(x, q.float(), a.eta)
            if a.criterion == "graph_geodesic":
                ck["graph"], cost = fit_graph(u, a.book, a.neighbors)
                unit = "squared_graph_length"
            else:
                cost = (
                    torch.cdist(
                        u, u[: a.book], compute_mode="donot_use_mm_for_euclid_dist"
                    )
                    .double()
                    .square()
                )
                unit = "squared_task_chord"
        lm, cert = solve(cost, a.beta, a.iterations, a.tolerance)
        cert["distortion_unit"] = unit
        if a.criterion == "posterior_kl":
            cert["distortion_bits"] = cert["distortion"]
        ck.update(
            kind="rd",
            alphabet=a.alphabet,
            outputs=b.cpu(),
            beta=a.beta,
            log_prior=lm,
            train_certificate=cert,
        )
        print(cert)
        save(a, ck)
    else:
        ck = torch.load(a.checkpoint, map_location="cpu", weights_only=True)
        if ck.get("kind") != "rd":
            raise ValueError("not an RD checkpoint")
        x, q, y, ids = test_source(a, ck, dev)
        q = q.double()
        b = ck["outputs"].to(dev)
        criterion = ck.get("criterion", "posterior_kl")
        if criterion == "posterior_kl":
            cost = kl_table(q, b)
        else:
            u = coordinates(x, q.float(), ck["eta"])
            if criterion == "graph_geodesic":
                cost = attached(u, ck["graph"])
            else:
                bx = ck["anchors"].to(dev)
                v = coordinates(bx, post(bx, ck), ck["eta"])
                cost = (
                    torch.cdist(u, v, compute_mode="donot_use_mm_for_euclid_dist")
                    .double()
                    .square()
                )
        W = (ck["log_prior"].to(dev) - LN2 * ck["beta"] * cost).softmax(1)
        report(
            a,
            ck,
            q,
            y,
            ids,
            W,
            b,
            dict(
                criterion=criterion,
                expected_coding_distortion=float((W * cost).sum(1).mean()),
                graph_k=ck["graph"]["k"] if criterion == "graph_geodesic" else None,
                beta=ck["beta"],
                alphabet=ck["alphabet"],
                train_certificate=ck["train_certificate"],
            ),
            cosine=x @ ck["anchors"].to(dev).T,
            fixed=cost if criterion != "posterior_kl" else None,
        )


if __name__ == "__main__":
    main()
