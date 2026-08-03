"""Tests for recipe loading and schema validation."""

import textwrap

import pytest

from migrator.config import load_recipe


def _write_recipe(tmp_path, body: str):
    recipe_dir = tmp_path / 'recipes'
    recipe_dir.mkdir()
    src_dir = tmp_path / 'src'
    src_dir.mkdir()
    recipe = recipe_dir / 'r.yaml'
    recipe.write_text(textwrap.dedent(body))
    return recipe


def test_minimal_recipe_loads(tmp_path):
    recipe = _write_recipe(tmp_path, """
        library: blas
        language: fortran
        source_dir: src
    """)
    config = load_recipe(recipe, project_root=tmp_path)
    assert config.library == 'blas'
    assert config.language == 'fortran'
    assert config.source_dir == tmp_path / 'src'


def test_missing_required_key_raises(tmp_path):
    recipe = _write_recipe(tmp_path, """
        library: blas
        # language and source_dir intentionally missing
    """)
    with pytest.raises(KeyError, match='language'):
        load_recipe(recipe, project_root=tmp_path)


def test_unknown_keys_warn_on_stderr(tmp_path, capsys):
    recipe = _write_recipe(tmp_path, """
        library: blas
        language: fortran
        source_dir: src
        copy-files: [FOO]   # typo: dash vs underscore
        spurious_key: 42
    """)
    config = load_recipe(recipe, project_root=tmp_path)
    captured = capsys.readouterr()
    assert 'unknown recipe key' in captured.err
    # The typo'd key should be mentioned by name, so the user can find
    # and fix the offending line.
    assert 'copy-files' in captured.err
    assert 'spurious_key' in captured.err
    # The mistyped key cannot have populated copy_files; that's the
    # exact silent-default-empty bug this warning guards against.
    assert config.copy_files == set()


# Both fields are retired — ``symbols`` (an nm-based archive scan) along
# with ``library_path``, ``prefix`` in v0.3.3 — but a recipe still carrying
# the key with an empty body must load without raising, so loader tolerance
# stays a regression-worthy property.
@pytest.mark.parametrize('key', ['symbols', 'prefix'])
def test_empty_retired_block_loads(tmp_path, key):
    recipe = _write_recipe(tmp_path, f"""
        library: blas
        language: fortran
        source_dir: src
        {key}:
    """)
    load_recipe(recipe, project_root=tmp_path)


def test_non_mapping_top_level_raises(tmp_path):
    recipe = _write_recipe(tmp_path, """
        - not_a_mapping
    """)
    with pytest.raises(ValueError, match='must be a mapping'):
        load_recipe(recipe, project_root=tmp_path)


@pytest.mark.parametrize('key', ['asymmetric_patches', 'one_sided_cleanup'])
def test_patch_classification_of_unlisted_patch_raises(tmp_path, key):
    # verify_patches() unions both classification lists into a skip set
    # without checking membership in patches:, so a name left behind
    # after its patch is dropped would silently weaken the CI check.
    recipe = _write_recipe(tmp_path, f"""
        library: blas
        language: fortran
        source_dir: src
        patches: [a.patch]
        {key}: [a.patch, gone.patch]
    """)
    with pytest.raises(ValueError, match='gone.patch'):
        load_recipe(recipe, project_root=tmp_path)


def test_skip_files_manifest_unions_with_inline_list(tmp_path):
    recipe = _write_recipe(tmp_path, """
        library: blas
        language: fortran
        source_dir: src
        skip_files: [dsdot]
        skip_files_manifest: skip.txt
    """)
    (recipe.parent / 'skip.txt').write_text(
        '# generated\n\nBLAS_cdot_c_s\nBLAS_cdot_s_c\n')
    config = load_recipe(recipe, project_root=tmp_path)
    assert config.skip_files == {'DSDOT', 'BLAS_CDOT_C_S', 'BLAS_CDOT_S_C'}


def _ep_recipe(tmp_path, symbols: str):
    recipe = _write_recipe(tmp_path, f"""
        library: mumps
        language: fortran
        source_dir: src
        ep_renames:
          manifest: privatize.txt
          symbols: {symbols}
    """)
    (recipe.parent / 'privatize.txt').write_text(
        '# exact linker-level names\ndescinit_\nblacs_gridinit_\n')
    return recipe


def test_ep_renames_resolve_against_manifest(tmp_path):
    config = load_recipe(_ep_recipe(tmp_path, '[DESCINIT, BLACS_GRIDINIT]'),
                         project_root=tmp_path)
    assert config.extra_renames == {'DESCINIT': 'EP_DESCINIT',
                                    'BLACS_GRIDINIT': 'EP_BLACS_GRIDINIT'}


def test_ep_rename_absent_from_manifest_raises(tmp_path):
    # The point of the key: the manifest owns the scheme, so a symbol it
    # does not privatize fails the load instead of surfacing as an
    # undefined ep_descset_ at final executable link.
    with pytest.raises(ValueError, match='DESCSET is not privatized'):
        load_recipe(_ep_recipe(tmp_path, '[DESCSET]'), project_root=tmp_path)
