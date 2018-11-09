// 2018 © William Chèvremont <william.chevremont@univ-grenoble-alpes.fr>

#pragma once

#include <boost/multi_array.hpp>
#include <pkg/dem/ScGeom.hpp>
#include <pkg/common/PeriodicEngines.hpp>


class PDFEngine : public PeriodicEngine {

public:
	class PDFCalculator {
	public: 
		PDFCalculator(string const& n) : name(n) {};
		virtual ~PDFCalculator() {};
		
		virtual vector<string> getSuffixes() const { return vector<string>().push_back(""); }
		virtual vector<string> getDatas() const = 0;
		virtual void cleanData() = 0;
		virtual bool addData(const shared_ptr<Interaction>&, Real const& dS ,Real const& V, int const& N) = 0;
		
		string name;
	};
	
	typedef boost::multi_array<shared_ptr<PDFCalculator>, 2> PDF;
	
	static void getSpectrums(vector<PDF> &);
	virtual void action() { LOG_WARN("Don't use this class directly! Use Law-specific implementation instead."); warnedOnce = true;};
	
	YADE_CLASS_BASE_DOC_ATTRS_CTOR_PY(PDFEngine, PeriodicEngine,
		"Base class for spectrums calculations. Should not be used directly. Compute Probability Density Functions in spherical coordinates and write result to a file.",
		((uint, numDiscretizeAngleTheta, 20,,"Number of sector for theta-angle"))
		((uint, numDiscretizeAnglePhi, 40,,"Number of sector for phi-angle"))
		//((Real, discretizeRadius, 0.1,,"d/a interval size"))
		((string, filename, "PDF.txt", , "Filename"))
		((bool, firstRun, true, (Attr::hidden | Attr::readonly), ""))
		((bool, warnedOnce, false, , "For one-time warning. May trigger usefull warnings"))
		,,
		//.def("getSpectrums", &LubricationDPFEngine::PyGetSpectrums,(py::arg("nPhi")=40, py::arg("nTheta")=20), "Get Stress spectrums").staticmethod("getSpectrums")
	);
	DECLARE_LOGGER;
	
	protected:
		void writeToFile(vector<PDF> const& );
};

REGISTER_SERIALIZABLE(PDFEngine);

template<class Phys>
class PDFSpheresStressCalculator : public PDFEngine::PDFCalculator {
public:
	PDFSpheresStressCalculator(Vector3r Phys::* member, string name) : PDFEngine::PDFCalculator(name), m_member(member), m_stress(Matrix3r::Zero()) {};
	vector<string> getSuffixes() const { return vector<string>({"xx", "xy", "xz", "yx", "yy", "yz", "zx", "zy", "zz"}); }
	vector<string> getDatas() const {
		vector<string> out;
		for(int i(0);i<3;i++) for(int j(0);j<3;j++) out.push_back(std::to_string(m_stress(i,j)));
		return out;
	}
	void cleanData() { m_stress = Matrix3r::Zero(); }
	bool addData(const shared_ptr<Interaction>& I, Real const& dS, Real const& V, int const& N) {
		if(!I->isReal()) return false;
	ScGeom* geom=YADE_CAST<ScGeom*>(I->geom.get());
	Phys* phys=YADE_CAST<Phys*>(I->phys.get());
	
	if(geom && phys) {
		Real r = geom->radius1 + geom->radius2 - geom->penetrationDepth;
		Vector3r l = r/(V*dS)*geom->normal;
		m_stress += phys->*(m_member)*l.transpose();
		return true;
	}
	else
		return false;
	}
	
private:
	Vector3r Phys::*m_member;
	Matrix3r m_stress;
};

class PDFSpheresDistanceCalculator : public PDFEngine::PDFCalculator {
	PDFSpheresDistanceCalculator(string name);
	vector<string> getDatas() const;
	void cleanData();
	bool addData(const shared_ptr<Interaction>&, Real const& dS ,Real const& V, int const& N);
	
private:
	Real m_h;
	uint m_N;
};

class PDFSpheresVelocityCalculator : public PDFEngine::PDFCalculator {
	PDFSpheresVelocityCalculator(string name);
	vector<string> getSuffixes() const;
	vector<string> getDatas() const;
	void cleanData();
	bool addData(const shared_ptr<Interaction>&, Real const& dS ,Real const& V, int const& N);
	
private:
	Vector3r m_vel;
	uint m_N;
};

class PDFSpheresIntrsCalculator : public PDFEngine::PDFCalculator {
	PDFSpheresIntrsCalculator(string name);
	vector<string> getDatas() const;
	void cleanData();
	bool addData(const shared_ptr<Interaction>&, Real const& dS ,Real const& V, int const& N);
	
private:
	Ream m_P;
};
