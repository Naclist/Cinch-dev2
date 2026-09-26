from cinch_dev2.workflow import load_workflow, verify_frozen_files


def test_workflow_has_all_registered_stages():
    stages = load_workflow()["stages"]
    assert len(stages) == 23  # numbered stages 01-22 plus the registered 10b substage
    assert len({stage["id"] for stage in stages}) == len(stages)
    assert stages[-1]["id"] == "22-diff-gwes"


def test_frozen_sources_are_unchanged():
    rows = verify_frozen_files()
    assert rows
    assert {row["status"] for row in rows} == {"OK"}
