"""Check exact rational interval uniqueness; not a period enclosure prover."""
import argparse
from fractions import Fraction
import json
from pathlib import Path
import sys


def rational(value):
    if not isinstance(value, str):
        raise ValueError("Endpoints and candidates must be strings, not floats")
    return Fraction(value)


def unique_integer(lower, upper):
    a, b = rational(lower), rational(upper)
    if a > b:
        raise ValueError("Reversed interval")
    first = -((-a.numerator) // a.denominator)
    last = b.numerator // b.denominator
    if first != last:
        raise ValueError("Interval does not contain exactly one integer")
    return first


def check(data):
    n = data.get("dimension")
    if type(n) is not int or n < 1:
        raise ValueError("Positive integer dimension required")
    candidate, intervals = data.get("candidate"), data.get("intervals")
    for matrix in (candidate, intervals):
        if not isinstance(matrix, list) or len(matrix) != n:
            raise ValueError("Matrix dimension mismatch")
        if any(not isinstance(row, list) or len(row) != n for row in matrix):
            raise ValueError("Matrix must be square")
    for i in range(n):
        for j in range(n):
            entry = intervals[i][j]
            if not isinstance(entry, list) or len(entry) != 2:
                raise ValueError(f"Invalid interval at {i+1},{j+1}")
            expected = rational(candidate[i][j])
            if expected.denominator != 1:
                raise ValueError("Candidate entry is not an integer")
            actual = unique_integer(*entry)
            if actual != expected:
                raise ValueError(f"Candidate mismatch at {i+1},{j+1}")
    return {"status": "interval_checks_passed", "entries": n*n,
            "scope": "Input enclosure validity, geometric integrality and saturation are not verified"}


def self_test():
    assert unique_integer("3.999", "4.001") == 4
    assert unique_integer("-4.001", "-3.999") == -4
    assert unique_integer("0", "0") == 0
    assert unique_integer("1/3", "1") == 1
    for pair in (("4.1", "4.2"), ("0", "1"), ("2", "1"), (0.1, "1")):
        try:
            unique_integer(*pair)
        except (ValueError, TypeError):
            pass
        else:
            raise AssertionError(f"Invalid interval accepted: {pair}")
    good = {"dimension": 1, "candidate": [["4"]], "intervals": [[["3.99", "4.01"]]]}
    assert check(good)["entries"] == 1
    for bad in (
        {**good, "candidate": [["5"]]},
        {**good, "candidate": [["1/2"]]},
        {**good, "intervals": []},
        {**good, "dimension": True},
    ):
        try:
            check(bad)
        except ValueError:
            pass
        else:
            raise AssertionError("Invalid matrix accepted")
    return {"status": "self_test_passed"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", nargs="?", type=Path)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    try:
        if args.self_test:
            result = self_test()
        elif args.input:
            result = check(json.loads(args.input.read_text(encoding="utf-8-sig")))
        else:
            parser.error("Provide an input JSON or --self-test")
        print(json.dumps(result, ensure_ascii=False))
        return 0
    except (ValueError, TypeError, OSError, KeyError) as error:
        print(json.dumps({"status": "failed", "reason": str(error)}, ensure_ascii=False))
        return 1


if __name__ == "__main__":
    sys.exit(main())
