@test convert_unit(10, orig="meV", to="GHz") == 10 * (1.0e-3 * 1.602176565e-19 / (1.0e9 * 6.62606957e-34))
@test qutip[:utilities][:convert_meV_to_GHz](10) == qt.utilities[:convert_meV_to_GHz](10)
