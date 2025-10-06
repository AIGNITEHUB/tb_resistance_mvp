"""
Quick tests for demo version
Run: python test_demo.py
"""

def test_imports():
    """Test all demo modules can be imported"""
    print("Testing imports...")
    try:
        from tb_resistance_mvp.drug_filtering_demo import (
            MutationGroupClassifier, 
            FakeMoleculeDatabase,
            FilterCriteriaHelper
        )
        print("✅ drug_filtering_demo imports OK")
    except ImportError as e:
        print(f"❌ Import failed: {e}")
        return False
    
    try:
        from tb_resistance_mvp.pipeline import mutational_scan
        print("✅ pipeline imports OK")
    except ImportError as e:
        print(f"❌ Import failed: {e}")
        return False
    
    return True


def test_classification():
    """Test mutation classification"""
    print("\nTesting mutation classification...")
    from tb_resistance_mvp.drug_filtering_demo import MutationGroupClassifier
    
    test_cases = [
        # (features, wt, mut, expected_group, description)
        (
            {"d_hydro": -3.5, "d_vol": 0.2, "d_charge": 0.0, "rel_pos": 0.45, "center_dist": 0.1, "motif_pen": 0.0},
            "R", "I", "H", "Hydrophobic increase (R→I)"
        ),
        (
            {"d_hydro": 5.0, "d_vol": -0.3, "d_charge": 0.0, "rel_pos": 0.35, "center_dist": 0.3, "motif_pen": 0.0},
            "I", "S", "P", "Polar increase (I→S)"
        ),
        (
            {"d_hydro": 0.0, "d_vol": 0.4, "d_charge": 2.0, "rel_pos": 0.5, "center_dist": 0.0, "motif_pen": 0.0},
            "D", "K", "C+", "Positive charge gain (D→K)"
        ),
        (
            {"d_hydro": 0.0, "d_vol": 0.0, "d_charge": -2.0, "rel_pos": 0.4, "center_dist": 0.2, "motif_pen": 0.0},
            "K", "E", "C-", "Negative charge gain (K→E)"
        ),
        (
            {"d_hydro": 1.0, "d_vol": -0.5, "d_charge": 0.0, "rel_pos": 0.4, "center_dist": 0.2, "motif_pen": 0.0},
            "F", "A", "A", "Aromatic loss (F→A)"
        ),
        (
            {"d_hydro": 0.0, "d_vol": 0.9, "d_charge": 0.0, "rel_pos": 0.3, "center_dist": 0.4, "motif_pen": 0.0},
            "G", "W", "V", "Volume increase (G→W)"
        ),
        (
            {"d_hydro": 0.5, "d_vol": 0.2, "d_charge": 0.0, "rel_pos": 0.5, "center_dist": 0.0, "motif_pen": 0.8},
            "L", "P", "M", "Motif penalty (L→P)"
        ),
    ]
    
    passed = 0
    for features, wt, mut, expected, desc in test_cases:
        result = MutationGroupClassifier.classify(features, wt, mut)
        if result == expected:
            print(f"  ✅ {desc}: {result}")
            passed += 1
        else:
            print(f"  ❌ {desc}: got {result}, expected {expected}")
    
    print(f"\nPassed {passed}/{len(test_cases)} classification tests")
    return passed == len(test_cases)


def test_molecule_search():
    """Test fake molecule database"""
    print("\nTesting molecule search...")
    from tb_resistance_mvp.drug_filtering_demo import FakeMoleculeDatabase
    
    groups = ["H", "P", "C+", "C-", "V", "A", "M"]
    
    for group in groups:
        molecules = FakeMoleculeDatabase.search(group, max_results=5)
        if len(molecules) > 0:
            print(f"  ✅ Group {group}: Found {len(molecules)} molecules")
            # Show first molecule
            mol = molecules[0]
            print(f"     Example: {mol['name']} (MW={mol['mw']:.1f}, LogP={mol['logp']:.2f})")
        else:
            print(f"  ❌ Group {group}: No molecules found")
            return False
    
    return True


def test_molecule_validation():
    """Test SMILES validation"""
    print("\nTesting molecule validation...")
    from tb_resistance_mvp.drug_filtering_demo import FakeMoleculeDatabase
    
    test_cases = [
        ("c1ccccc1", "H", "Benzene for hydrophobic"),
        ("CCO", "P", "Ethanol for polar"),
        ("c1ccc(cc1)C(=O)O", "C+", "Benzoic acid for C+"),
    ]
    
    for smiles, group, desc in test_cases:
        result = FakeMoleculeDatabase.validate_molecule(smiles, group)
        status = "✅" if result["pass"] else "⚠️"
        print(f"  {status} {desc}")
        print(f"     MW={result['scores']['MW']:.1f}, LogP={result['scores']['LogP']:.2f}")
        if result['violations']:
            print(f"     Violations: {', '.join(result['violations'])}")
    
    return True


def test_filter_descriptions():
    """Test filter criteria descriptions"""
    print("\nTesting filter descriptions...")
    from tb_resistance_mvp.drug_filtering_demo import FilterCriteriaHelper
    
    groups = ["H", "P", "C+", "C-", "V", "A", "M"]
    
    for group in groups:
        desc = FilterCriteriaHelper.get_description(group)
        if len(desc) > 50:  # Should have substantial content
            print(f"  ✅ Group {group}: {desc.split(chr(10))[0][:60]}...")
        else:
            print(f"  ❌ Group {group}: Description too short")
            return False
    
    return True


def test_full_pipeline():
    """Test complete mutational scan with demo mode"""
    print("\nTesting full pipeline...")
    from tb_resistance_mvp.pipeline import mutational_scan
    
    sequence = "MASTKQLLAVGHVPRNTLDEYQFGL" * 10  # 250 AA
    
    try:
        result = mutational_scan(
            sequence=sequence,
            region_start=50,
            region_end=70,
            gene_hint="rpoB",
            enable_drug_filtering=True,
            resistance_cutoff=0.7,
            demo_mode=True
        )
        
        # Check structure
        assert "matrix" in result, "Missing matrix"
        assert "drug_filtering" in result, "Missing drug_filtering"
        assert "heatmap" in result, "Missing heatmap"
        
        # Check drug filtering
        df = result["drug_filtering"]
        assert df["enabled"] == True, "Drug filtering not enabled"
        assert df["demo_mode"] == True, "Demo mode not enabled"
        assert "filter_summary" in df, "Missing filter_summary"
        assert "example_molecules" in df, "Missing example_molecules"
        
        # Check mutations have groups
        for mut in result["matrix"][:10]:  # Check first 10
            assert "group" in mut, f"Missing group in mutation {mut}"
            assert mut["group"] in ["H", "P", "C+", "C-", "V", "A", "M", "General"], \
                f"Invalid group: {mut['group']}"
        
        print(f"  ✅ Pipeline executed successfully")
        print(f"     - Scanned {len(result['matrix'])} mutations")
        print(f"     - Found {df['high_risk_count']} high-risk mutations")
        print(f"     - Identified {len(df['filter_summary'])} groups")
        
        # Show group distribution
        from collections import Counter
        groups = [m["group"] for m in result["matrix"] if max(m["p_resist"].values()) >= 0.7]
        group_counts = Counter(groups)
        print(f"     - Group distribution: {dict(group_counts)}")
        
        return True
        
    except Exception as e:
        print(f"  ❌ Pipeline failed: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_example_molecules():
    """Test that example molecules are returned for each group"""
    print("\nTesting example molecules...")
    from tb_resistance_mvp.pipeline import mutational_scan
    
    sequence = "MASTKQLLAVGHVPRNTLDEYQFGL" * 10
    
    result = mutational_scan(
        sequence=sequence,
        region_start=100,
        region_end=120,
        gene_hint="rpoB",
        enable_drug_filtering=True,
        resistance_cutoff=0.6,  # Lower to get more groups
        demo_mode=True
    )
    
    if "example_molecules" in result["drug_filtering"]:
        example_mols = result["drug_filtering"]["example_molecules"]
        
        for group, molecules in example_mols.items():
            if len(molecules) > 0:
                print(f"  ✅ Group {group}: {len(molecules)} example molecules")
                # Show first molecule
                mol = molecules[0]
                print(f"     {mol['name']}: {mol['smiles']}")
            else:
                print(f"  ⚠️ Group {group}: No molecules (group might have no high-risk mutations)")
        
        return len(example_mols) > 0
    else:
        print("  ❌ No example molecules found")
        return False


def run_all_tests():
    """Run all tests"""
    print("=" * 60)
    print("TB Resistance Predictor - Demo Version Tests")
    print("=" * 60)
    
    tests = [
        ("Imports", test_imports),
        ("Classification", test_classification),
        ("Molecule Search", test_molecule_search),
        ("Molecule Validation", test_molecule_validation),
        ("Filter Descriptions", test_filter_descriptions),
        ("Full Pipeline", test_full_pipeline),
        ("Example Molecules", test_example_molecules),
    ]
    
    results = {}
    for name, test_func in tests:
        print(f"\n{'='*60}")
        try:
            results[name] = test_func()
        except Exception as e:
            print(f"❌ {name} raised exception: {e}")
            import traceback
            traceback.print_exc()
            results[name] = False
    
    # Summary
    print("\n" + "=" * 60)
    print("TEST SUMMARY")
    print("=" * 60)
    
    for name, passed in results.items():
        status = "✅ PASS" if passed else "❌ FAIL"
        print(f"{status:12} {name}")
    
    total_passed = sum(results.values())
    total_tests = len(results)
    
    print("=" * 60)
    print(f"Total: {total_passed}/{total_tests} tests passed")
    
    if total_passed == total_tests:
        print("\n🎉 All tests passed! Demo is ready to run.")
        print("\nNext steps:")
        print("  1. Run: streamlit run ui_streamlit_demo.py")
        print("  2. Open: http://localhost:8501")
        print("  3. Try scanning mutations with drug filtering enabled")
    else:
        print(f"\n⚠️ {total_tests - total_passed} test(s) failed. Check errors above.")
    
    return total_passed == total_tests


if __name__ == "__main__":
    import sys
    success = run_all_tests()
    sys.exit(0 if success else 1)