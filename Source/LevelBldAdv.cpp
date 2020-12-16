
#include <AMReX_LevelBld.H>
#include <projamrex.H>

using namespace amrex;

class LevelBldAdv
    :
    public LevelBld
{
    virtual void variableSetUp () override;
    virtual void variableCleanUp () override;
    virtual AmrLevel *operator() () override;
    virtual AmrLevel *operator() (Amr&            papa,
                                  int             lev,
                                  const Geometry& level_geom,
                                  const BoxArray& ba,
				  const DistributionMapping& dm,
                                  Real            time) override;
};

LevelBldAdv Adv_bld;

LevelBld*
getLevelBld ()
{
    return &Adv_bld;
}

void
LevelBldAdv::variableSetUp ()
{
    Projamrex::variableSetUp();
}

void
LevelBldAdv::variableCleanUp ()
{
    Projamrex::variableCleanUp();
}

AmrLevel*
LevelBldAdv::operator() ()
{
    return new Projamrex;
}

AmrLevel*
LevelBldAdv::operator() (Amr&            papa,
	   	         int             lev,
                         const Geometry& level_geom,
                         const BoxArray& ba,
                         const DistributionMapping& dm,
                         Real            time)
{
    return new Projamrex(papa, lev, level_geom, ba, dm, time);
}
