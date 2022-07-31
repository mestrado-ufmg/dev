from numpy import array, linspace, zeros_like
from pybird.models.enums import TailShape

from pybird.models.geo_model import GeoModel
from pybird.geo.utils.system import BaseSystem, TailSystem, WingSystem
from pybird.geo.utils import points, sections, surface, bezier
from pybird.models.types import Curve, Point

class Geo:

    def __init__(self, data: GeoModel) -> None:
        self.data = data
        self._baseSystem = BaseSystem()
        return

    def build(self) -> None:
        print('- Building geometry')
        self._build_wing()
        self._build_body()
        self._build_head()
        self._build_tail()
        return
    
    def _symmetrical_point(self, a: Point) -> Point:
        return array([a[0], -a[1], a[2]])
    
    def _symmetrical_curve(self, a: Curve) -> Curve:
        out = zeros_like(a)
        for i in range(len(a[:, 0])):
            out[i, :] = self._symmetrical_point(a[i, :])
        return out
    
    def _build_wing(self) -> None:
        
        # Left wing

        # Wing systems
        leftSystem = WingSystem(
            thetaRootY=self.data.wing.thetaRootY,
            theta1=self.data.wing.theta1_e,
            theta2=self.data.wing.theta2_e,
            theta3=self.data.wing.theta3_e,
            theta4=self.data.wing.theta4_e,
            theta5=self.data.wing.theta5_e,
            theta6=self.data.wing.theta6_e,
            theta7=self.data.wing.theta7_e,
            x0=self._baseSystem.x0,
            y0=self._baseSystem.y0,
            z0=self._baseSystem.z0,
        )

        # Countor points
        self.p0e = points.p0(self.data.wing, leftSystem)
        l1 = points.l1(self.data.wing, leftSystem, self.p0e)
        l2 = points.l2(self.data.wing, leftSystem, l1)
        l3 = points.l3(self.data.wing, leftSystem, l2)
        self.p1e = points.p1(self.data.wing, leftSystem, self.p0e, True)
        p3eLine = points.p3Line(self.data.wing, leftSystem, l2)
        self.aux1e = points.aux1(self.data.wing, leftSystem, l1, self.p1e, p3eLine)
        self.p2e = points.p2(self.p1e, self.aux1e, p3eLine)
        self.p3e = points.p3(self.data.wing, leftSystem, l2)
        self.p4e = points.p4(self.data.wing, leftSystem, l2)
        self.p5e = points.p5(self.data.wing, leftSystem, l3)
        self.p7e = points.p7(self.data.wing, leftSystem, self.p5e)
        tp6e = 0.8
        self.aux2e = points.aux2(self.data.wing, leftSystem, self.p5e, self.p7e)
        self.p6e = points.p6(self.p5e, self.aux2e, self.p7e, tp6e)
        p10_11_aux = points.p10_11_aux(self.data.wing, leftSystem, self.p4e, self.p3e)
        self.aux3e = points.aux3(self.data.wing, leftSystem, self.p7e, p10_11_aux)
        self.p10e = points.p10(self.data.wing, p10_11_aux, self.p7e, self.aux3e)
        tp8e = 0.2
        self.p8e = points.p8(self.p7e, self.aux3e, p10_11_aux, tp8e)
        tp9e = 0.6
        self.p9e = points.p9(self.p7e, self.aux3e, p10_11_aux, tp9e)
        self.p13e = points.p13(self.data.wing, leftSystem, self.p0e, True)
        self.aux4e = points.aux4(self.data.wing, leftSystem, p10_11_aux, self.p13e)
        self.p11e = points.p11(self.data.wing, p10_11_aux, self.p13e, self.aux4e)
        self.p12e = points.p12(p10_11_aux, self.aux4e, self.p13e)

        # Control points
        self.c1e = points.c1(self.p1e, self.aux1e, p3eLine)
        self.c2e = points.c2(self.data.wing, self.p1e, self.p2e, self.p3e)
        self.c3e = points.c3(self.data.wing, leftSystem, self.p2e, self.p3e, self.p4e)
        self.c4e = points.c4(self.data.wing, leftSystem, self.p3e, self.p4e)
        self.c5e = points.c5(self.data.wing, leftSystem, self.p3e, self.p4e)
        self.c6e, self.c7e = points.c6_c7(self.p5e, self.aux2e, self.p7e, tp6e)
        self.c8e, self.c9e = points.c8_c9(self.p5e, self.aux2e, self.p7e, tp6e)
        self.c10e, self.c11e = points.c10_c11(self.p7e, self.aux3e, self.p10e, tp8e)
        self.c12e, self.c13e = points.c12_c13(self.p7e, self.aux3e, self.p10e, tp8e, tp9e)
        self.c14e, self.c15e = points.c14_c15(self.p7e, self.aux3e, self.p10e, tp9e)
        self.c16e = points.c16(self.data.wing, self.p10e, self.p11e, self.c15e)
        self.c18e, self.c19e = points.c18_c19(self.p11e, self.aux4e, self.p13e, 0.5)
        self.c17e = points.c17(self.data.wing, self.p10e, self.p11e, self.c18e)
        self.c20e, self.c21e = points.c20_c21(self.p11e, self.aux4e, self.p13e, 0.5)

        # Sections
        self.p14e, self.p15e, self.curve13e, self.curve14e, self.curve15e, self.curve16e = sections.process_section(self.data.wing.foils[0], self.p1e, self.p13e, leftSystem.x1, leftSystem.z1, n=100)

        p16e1, p17e1, curve19e1, curve20e1, curve21e1, curve22e1 = sections.process_section(self.data.wing.foils[0], self.p2e, self.p12e, leftSystem.x1, leftSystem.z1, n=100)
        p16e2, p17e2, curve19e2, curve20e2, curve21e2, curve22e2 = sections.process_section(self.data.wing.foils[1], self.p2e, self.p12e, leftSystem.x1, leftSystem.z1, n=100)
        self.p16e, self.p17e = 0.5 * (p16e1 + p16e2), 0.5 * (p17e1 + p17e2)
        self.curve19e, self.curve20e = 0.5 * (curve19e1 + curve19e2), 0.5 * (curve20e1 + curve20e2)
        self.curve21e, self.curve22e = 0.5 * (curve21e1 + curve21e2), 0.5 * (curve22e1 + curve22e2)

        self.p18e, self.p19e, self.curve25e, self.curve26e, self.curve27e, self.curve28e = sections.process_section(self.data.wing.foils[1], self.p3e, self.p11e, leftSystem.x2Tip, leftSystem.z2Tip, n=100)
        self.p20e, self.p21e, self.curve31e, self.curve32e, self.curve33e, self.curve34e = sections.process_section(self.data.wing.foils[1], self.p4e, self.p10e, leftSystem.x3, leftSystem.z3, n=100)
        
        dataAux1 = sections.process_section(self.data.wing.foils[1], self.p5e, self.p9e, leftSystem.x3, leftSystem.z3, n=100)
        dataAux2 = sections.process_section(self.data.wing.foils[2], self.p5e, self.p9e, leftSystem.x3, leftSystem.z3, n=100)
        self.p22e, self.p23e, self.curve37e, self.curve38e, self.curve39e, self.curve40e = tp9e * dataAux1[0] + (1 - tp9e) * dataAux2[0], tp9e * dataAux1[1] + (1 - tp9e) * dataAux2[1], tp9e * dataAux1[2] + (1 - tp9e) * dataAux2[2], tp9e * dataAux1[3] + (1 - tp9e) * dataAux2[3], tp9e * dataAux1[4] + (1 - tp9e) * dataAux2[4], tp9e * dataAux1[5] + (1 - tp9e) * dataAux2[5]
        
        self.p24e, self.p25e, self.curve43e, self.curve44e, self.curve45e, self.curve46e = sections.process_section(self.data.wing.foils[1], self.p6e, self.p8e, leftSystem.x3, leftSystem.z3, n=100)

        # Curves between sections
        self.curve17e = surface.interpolate_curve(
            self.p14e,
            self.p16e,
            bezier.quadratic([self.p1e, self.c1e, self.p2e], linspace(0, 1, num=50)),
            bezier.cubic([self.p13e, self.c21e, self.c20e, self.p12e], linspace(0, 1, num=50)),
        )
        self.curve18e = surface.interpolate_curve(
            self.p15e,
            self.p17e,
            bezier.quadratic([self.p1e, self.c1e, self.p2e], linspace(0, 1, num=50)),
            bezier.cubic([self.p13e, self.c21e, self.c20e, self.p12e], linspace(0, 1, num=50)),
        )

        self.curve23e = surface.interpolate_curve(
            self.p16e,
            self.p18e,
            bezier.cubic([self.p2e, self.c2e, self.c3e, self.p3e], linspace(0, 1, num=50)),
            bezier.cubic([self.p12e, self.c19e, self.c18e, self.p11e], linspace(0, 1, num=50)),
        )
        self.curve24e = surface.interpolate_curve(
            self.p17e,
            self.p19e,
            bezier.cubic([self.p2e, self.c2e, self.c3e, self.p3e], linspace(0, 1, num=50)),
            bezier.cubic([self.p12e, self.c19e, self.c18e, self.p11e], linspace(0, 1, num=50)),
        )

        self.curve29e = surface.interpolate_curve(
            self.p18e,
            self.p20e,
            bezier.cubic([self.p3e, self.c4e, self.c5e, self.p4e], linspace(0, 1, num=50)),
            bezier.cubic([self.p11e, self.c17e, self.c16e, self.p10e], linspace(0, 1, num=50)),
        )
        self.curve30e = surface.interpolate_curve(
            self.p19e,
            self.p21e,
            bezier.cubic([self.p3e, self.c4e, self.c5e, self.p4e], linspace(0, 1, num=50)),
            bezier.cubic([self.p11e, self.c17e, self.c16e, self.p10e], linspace(0, 1, num=50)),
        )

        self.curve35e = surface.interpolate_curve(
            self.p20e,
            self.p22e,
            bezier.quadratic([self.p4e, 0.5 * (self.p4e + self.p5e), self.p5e], linspace(0, 1, num=50)),
            bezier.cubic([self.p10e, self.c15e, self.c14e, self.p9e], linspace(0, 1, num=50)),
        )
        self.curve36e = surface.interpolate_curve(
            self.p21e,
            self.p23e,
            bezier.quadratic([self.p4e, 0.5 * (self.p4e + self.p5e), self.p5e], linspace(0, 1, num=50)),
            bezier.cubic([self.p10e, self.c15e, self.c14e, self.p9e], linspace(0, 1, num=50)),
        )

        self.curve41e = surface.interpolate_curve(
            self.p22e,
            self.p24e,
            bezier.cubic([self.p5e, self.c6e, self.c7e, self.p6e], linspace(0, 1, num=50)),
            bezier.cubic([self.p9e, self.c13e, self.c12e, self.p8e], linspace(0, 1, num=50)),
        )
        self.curve42e = surface.interpolate_curve(
            self.p23e,
            self.p25e,
            bezier.cubic([self.p5e, self.c6e, self.c7e, self.p6e], linspace(0, 1, num=50)),
            bezier.cubic([self.p9e, self.c13e, self.c12e, self.p8e], linspace(0, 1, num=50)),
        )
        self.curve47e = surface.interpolate_tip_curve(self.p24e, self.p7e)
        self.curve48e = surface.interpolate_tip_curve(self.p25e, self.p7e)

        # Right wing
        
        # Wing systems
        rightSystem = WingSystem(
            thetaRootY=self.data.wing.thetaRootY,
            theta1=self.data.wing.theta1_d,
            theta2=self.data.wing.theta2_d,
            theta3=self.data.wing.theta3_d,
            theta4=self.data.wing.theta4_d,
            theta5=self.data.wing.theta5_d,
            theta6=self.data.wing.theta6_d,
            theta7=self.data.wing.theta7_d,
            x0=self._baseSystem.x0,
            y0=self._baseSystem.y0,
            z0=self._baseSystem.z0,
        )

        # Countor points
        p0d = points.p0(self.data.wing, rightSystem)
        l1 = points.l1(self.data.wing, rightSystem, p0d)
        l2 = points.l2(self.data.wing, rightSystem, l1)
        l3 = points.l3(self.data.wing, rightSystem, l2)
        p1d = points.p1(self.data.wing, rightSystem, p0d, True)
        p3dLine = points.p3Line(self.data.wing, rightSystem, l2)
        aux1d = points.aux1(self.data.wing, rightSystem, l1, p1d, p3dLine)
        p2d = points.p2(p1d, aux1d, p3dLine)
        p3d = points.p3(self.data.wing, rightSystem, l2)
        p4d = points.p4(self.data.wing, rightSystem, l2)
        p5d = points.p5(self.data.wing, rightSystem, l3)
        p7d = points.p7(self.data.wing, rightSystem, p5d)
        tp6d = 0.8
        aux2d = points.aux2(self.data.wing, rightSystem, p5d, p7d)
        p6d = points.p6(p5d, aux2d, p7d, tp6d)
        p10_11_aux = points.p10_11_aux(self.data.wing, rightSystem, p4d, p3d)
        aux3d = points.aux3(self.data.wing, rightSystem, p7d, p10_11_aux)
        p10d = points.p10(self.data.wing, p10_11_aux, p7d, aux3d)
        tp8d = 0.2
        p8d = points.p8(p7d, aux3d, p10_11_aux, tp8d)
        tp9d = 0.6
        p9d = points.p9(p7d, aux3d, p10_11_aux, tp9d)
        p13d = points.p13(self.data.wing, rightSystem, p0d, True)
        aux4d = points.aux4(self.data.wing, rightSystem, p10_11_aux, p13d)
        p11d = points.p11(self.data.wing, p10_11_aux, p13d, aux4d)
        p12d = points.p12(p10_11_aux, aux4d, p13d)

        # Control points
        c1d = points.c1(p1d, aux1d, p3dLine)
        c2d = points.c2(self.data.wing, p1d, p2d, p3d)
        c3d = points.c3(self.data.wing, rightSystem, p2d, p3d, p4d)
        c4d = points.c4(self.data.wing, rightSystem, p3d, p4d)
        c5d = points.c5(self.data.wing, rightSystem, p3d, p4d)
        c6d, c7d = points.c6_c7(p5d, aux2d, p7d, tp6d)
        c8d, c9d = points.c8_c9(p5d, aux2d, p7d, tp6d)
        c10d, c11d = points.c10_c11(p7d, aux3d, p10d, tp8d)
        c12d, c13d = points.c12_c13(p7d, aux3d, p10d, tp8d, tp9d)
        c14d, c15d = points.c14_c15(p7d, aux3d, p10d, tp9d)
        c16d = points.c16(self.data.wing, p10d, p11d, c15d)
        c18d, c19d = points.c18_c19(p11d, aux4d, p13d, 0.5)
        c17d = points.c17(self.data.wing, p10d, p11d, c18d)
        c20d, c21d = points.c20_c21(p11d, aux4d, p13d, 0.5)

        # Sections
        p14d, p15d, curve13d, curve14d, curve15d, curve16d = sections.process_section(self.data.wing.foils[0], p1d, p13d, -rightSystem.x1, rightSystem.z1, n=100)

        p16d1, p17d1, curve19d1, curve20d1, curve21d1, curve22d1 = sections.process_section(self.data.wing.foils[0], p2d, p12d, -rightSystem.x1, rightSystem.z1, n=100)
        p16d2, p17d2, curve19d2, curve20d2, curve21d2, curve22d2 = sections.process_section(self.data.wing.foils[1], p2d, p12d, -rightSystem.x1, rightSystem.z1, n=100)
        p16d, p17d = 0.5 * (p16d1 + p16d2), 0.5 * (p17d1 + p17d2)
        curve19d, curve20d = 0.5 * (curve19d1 + curve19d2), 0.5 * (curve20d1 + curve20d2)
        curve21d, curve22d = 0.5 * (curve21d1 + curve21d2), 0.5 * (curve22d1 + curve22d2)

        p18d, p19d, curve25d, curve26d, curve27d, curve28d = sections.process_section(self.data.wing.foils[1], p3d, p11d, rightSystem.x2Tip, rightSystem.z2Tip, n=100)
        p20d, p21d, curve31d, curve32d, curve33d, curve34d = sections.process_section(self.data.wing.foils[1], p4d, p10d, rightSystem.x3, rightSystem.z3, n=100)
        
        dataAux1 = sections.process_section(self.data.wing.foils[1], p5d, p9d, -rightSystem.x3, rightSystem.z3, n=100)
        dataAux2 = sections.process_section(self.data.wing.foils[2], p5d, p9d, -rightSystem.x3, rightSystem.z3, n=100)
        p22d, p23d, curve37d, curve38d, curve39d, curve40d = tp9d * dataAux1[0] + (1 - tp9d) * dataAux2[0], tp9d * dataAux1[1] + (1 - tp9d) * dataAux2[1], tp9d * dataAux1[2] + (1 - tp9d) * dataAux2[2], tp9d * dataAux1[3] + (1 - tp9d) * dataAux2[3], tp9d * dataAux1[4] + (1 - tp9d) * dataAux2[4], tp9d * dataAux1[5] + (1 - tp9d) * dataAux2[5]
        
        p24d, p25d, curve43d, curve44d, curve45d, curve46d = sections.process_section(self.data.wing.foils[1], p6d, p8d, -rightSystem.x3, rightSystem.z3, n=100)

        # Curves between sections
        curve17d = surface.interpolate_curve(
            p14d,
            p16d,
            bezier.quadratic([p1d, c1d, p2d], linspace(0, 1, num=50)),
            bezier.cubic([p13d, c21d, c20d, p12d], linspace(0, 1, num=50)),
        )
        curve18d = surface.interpolate_curve(
            p15d,
            p17d,
            bezier.quadratic([p1d, c1d, p2d], linspace(0, 1, num=50)),
            bezier.cubic([p13d, c21d, c20d, p12d], linspace(0, 1, num=50)),
        )

        curve23d = surface.interpolate_curve(
            p16d,
            p18d,
            bezier.cubic([p2d, c2d, c3d, p3d], linspace(0, 1, num=50)),
            bezier.cubic([p12d, c19d, c18d, p11d], linspace(0, 1, num=50)),
        )
        curve24d = surface.interpolate_curve(
            p17d,
            p19d,
            bezier.cubic([p2d, c2d, c3d, p3d], linspace(0, 1, num=50)),
            bezier.cubic([p12d, c19d, c18d, p11d], linspace(0, 1, num=50)),
        )

        curve29d = surface.interpolate_curve(
            p18d,
            p20d,
            bezier.cubic([p3d, c4d, c5d, p4d], linspace(0, 1, num=50)),
            bezier.cubic([p11d, c17d, c16d, p10d], linspace(0, 1, num=50)),
        )
        curve30d = surface.interpolate_curve(
            p19d,
            p21d,
            bezier.cubic([p3d, c4d, c5d, p4d], linspace(0, 1, num=50)),
            bezier.cubic([p11d, c17d, c16d, p10d], linspace(0, 1, num=50)),
        )

        curve35d = surface.interpolate_curve(
            p20d,
            p22d,
            bezier.quadratic([p4d, 0.5 * (p4d + p5d), p5d], linspace(0, 1, num=50)),
            bezier.cubic([p10d, c15d, c14d, p9d], linspace(0, 1, num=50)),
        )
        curve36d = surface.interpolate_curve(
            p21d,
            p23d,
            bezier.quadratic([p4d, 0.5 * (p4d + p5d), p5d], linspace(0, 1, num=50)),
            bezier.cubic([p10d, c15d, c14d, p9d], linspace(0, 1, num=50)),
        )

        curve41d = surface.interpolate_curve(
            p22d,
            p24d,
            bezier.cubic([p5d, c6d, c7d, p6d], linspace(0, 1, num=50)),
            bezier.cubic([p9d, c13d, c12d, p8d], linspace(0, 1, num=50)),
        )
        curve42d = surface.interpolate_curve(
            p23d,
            p25d,
            bezier.cubic([p5d, c6d, c7d, p6d], linspace(0, 1, num=50)),
            bezier.cubic([p9d, c13d, c12d, p8d], linspace(0, 1, num=50)),
        )
        curve47d = surface.interpolate_tip_curve(p24d, p7d)
        curve48d = surface.interpolate_tip_curve(p25d, p7d)

        # Symetrical points and curves
        self.p0d = self._symmetrical_point(p0d)
        self.p1d = self._symmetrical_point(p1d)
        self.aux1d = self._symmetrical_point(aux1d)
        self.p2d = self._symmetrical_point(p2d)
        self.p3d = self._symmetrical_point(p3d)
        self.p4d = self._symmetrical_point(p4d)
        self.p5d = self._symmetrical_point(p5d)
        self.p7d = self._symmetrical_point(p7d)
        self.aux2d = self._symmetrical_point(aux2d)
        self.p6d = self._symmetrical_point(p6d)
        self.aux3d = self._symmetrical_point(aux3d)
        self.p10d = self._symmetrical_point(p10d)
        self.p8d = self._symmetrical_point(p8d)
        self.p9d = self._symmetrical_point(p9d)
        self.p13d = self._symmetrical_point(p13d)
        self.aux4d = self._symmetrical_point(aux4d)
        self.p11d = self._symmetrical_point(p11d)
        self.p12d = self._symmetrical_point(p12d)

        self.c1d = self._symmetrical_point(c1d)
        self.c2d = self._symmetrical_point(c2d)
        self.c3d = self._symmetrical_point(c3d)
        self.c4d = self._symmetrical_point(c4d)
        self.c5d = self._symmetrical_point(c5d)
        self.c6d, self.c7d = self._symmetrical_point(c6d), self._symmetrical_point(c7d)
        self.c8d, self.c9d = self._symmetrical_point(c8d), self._symmetrical_point(c9d)
        self.c10d, self.c11d = self._symmetrical_point(c10d), self._symmetrical_point(c11d)
        self.c12d, self.c13d = self._symmetrical_point(c12d), self._symmetrical_point(c13d)
        self.c14d, self.c15d = self._symmetrical_point(c14d), self._symmetrical_point(c15d)
        self.c16d = self._symmetrical_point(c16d)
        self.c18d, self.c19d = self._symmetrical_point(c18d), self._symmetrical_point(c19d)
        self.c17d = self._symmetrical_point(c17d)
        self.c20d, self.c21d = self._symmetrical_point(c20d), self._symmetrical_point(c21d)

        self.p14d, self.p15d = self._symmetrical_point(p14d), self._symmetrical_point(p15d)
        self.curve13d, self.curve14d, self.curve15d, self.curve16d = self._symmetrical_curve(curve13d), self._symmetrical_curve(curve14d), self._symmetrical_curve(curve15d), self._symmetrical_curve(curve16d)

        self.p16d, self.p17d = self._symmetrical_point(p16d), self._symmetrical_point(p17d)
        self.curve19d, self.curve20d, self.curve21d, self.curve22d = self._symmetrical_curve(curve19d), self._symmetrical_curve(curve20d), self._symmetrical_curve(curve21d), self._symmetrical_curve(curve22d)

        self.p18d, self.p19d = self._symmetrical_point(p18d), self._symmetrical_point(p19d)
        self.curve25d, self.curve26d, self.curve27d, self.curve28d = self._symmetrical_curve(curve25d), self._symmetrical_curve(curve26d), self._symmetrical_curve(curve27d), self._symmetrical_curve(curve28d)

        self.p20d, self.p21d = self._symmetrical_point(p20d), self._symmetrical_point(p21d)
        self.curve31d, self.curve32d, self.curve33d, self.curve34d = self._symmetrical_curve(curve31d), self._symmetrical_curve(curve32d), self._symmetrical_curve(curve33d), self._symmetrical_curve(curve34d)

        self.p22d, self.p23d = self._symmetrical_point(p22d), self._symmetrical_point(p23d)
        self.curve37d, self.curve38d, self.curve39d, self.curve40d = self._symmetrical_curve(curve37d), self._symmetrical_curve(curve38d), self._symmetrical_curve(curve39d), self._symmetrical_curve(curve40d)

        self.p24d, self.p25d = self._symmetrical_point(p24d), self._symmetrical_point(p25d)
        self.curve43d, self.curve44d, self.curve45d, self.curve46d = self._symmetrical_curve(curve43d), self._symmetrical_curve(curve44d), self._symmetrical_curve(curve45d), self._symmetrical_curve(curve46d)

        self.curve17d = self._symmetrical_curve(curve17d)
        self.curve18d = self._symmetrical_curve(curve18d)
        self.curve23d = self._symmetrical_curve(curve23d)
        self.curve24d = self._symmetrical_curve(curve24d)
        self.curve29d = self._symmetrical_curve(curve29d)
        self.curve30d = self._symmetrical_curve(curve30d)
        self.curve35d = self._symmetrical_curve(curve35d)
        self.curve36d = self._symmetrical_curve(curve36d)
        self.curve41d = self._symmetrical_curve(curve41d)
        self.curve42d = self._symmetrical_curve(curve42d)
        self.curve47d = self._symmetrical_curve(curve47d)
        self.curve48d = self._symmetrical_curve(curve48d)
        
        return
    
    def _build_body(self) -> None:

        self.p26e = points.p26(self.data.wing, self.data.body, self._baseSystem)
        self.p27e = points.p27(self.data.wing, self.data.body, self._baseSystem)
        self.p26d = self._symmetrical_point(self.p26e)
        self.p27d = self._symmetrical_point(self.p27e)
        self.p28 = points.p28(self.data.wing, self.data.body, self._baseSystem)
        self.p29 = points.p29(self.data.wing, self.data.body, self._baseSystem)
        self.p31 = points.p31(self.data.wing, self.data.body, self._baseSystem)
        aux5 = points.aux5(self.data.body, self._baseSystem, self.p29, self.p31)
        aux6 = points.aux6(self.data.body, self._baseSystem, self.p29, self.p31)
        self.p30, t30 = points.p30(self.p29, aux5, aux6, self.p31, self.p14e)
        self.p32 = points.p32(self.data.wing, self.data.body, self._baseSystem)
        self.p33 = points.p33(self.data.wing, self.data.body, self._baseSystem)
        self.p34 = points.p34(self.data.wing, self.data.body, self._baseSystem)
        self.p36 = points.p36(self.data.wing, self.data.body, self._baseSystem)
        aux7 = points.aux7(self.data.body, self._baseSystem, self.p34, self.p36)
        aux8 = points.aux8(self.data.body, self._baseSystem, self.p34, self.p36)
        self.p35, t35 = points.p35(self.p34, aux7, aux8, self.p36, self.p14e)
        self.p37 = points.p37(self.data.wing, self.data.body, self._baseSystem)

        self.c22e = points.c22(self.data.body, self._baseSystem, self.p26e, self.p1e)
        self.c23e = points.c23(self.data.body, self._baseSystem, self.p26e, self.p1e)
        self.c24e = points.c24(self.data.body, self._baseSystem, self.p13e, self.p27e)
        self.c25e = points.c25(self.data.body, self._baseSystem, self.p13e, self.p27e)
        self.c22d = self._symmetrical_point(self.c22e)
        self.c23d = self._symmetrical_point(self.c23e)
        self.c24d = self._symmetrical_point(self.c24e)
        self.c25d = self._symmetrical_point(self.c25e)
        self.c26 = points.c26(self.data.body, self._baseSystem, self.p28, self.p29)
        self.c27 = points.c27(self.data.body, self._baseSystem, self.p28, self.p29)
        self.c28, self.c29 = points.c28_29(self.p29, aux5, aux6, self.p31, t30)
        self.c30, self.c31 = points.c30_31(self.p29, aux5, aux6, self.p31, t30)
        self.c32 = points.c32(self.data.body, self._baseSystem, self.p31, self.p32)
        self.c33 = points.c33(self.data.body, self._baseSystem, self.p31, self.p32)
        self.c34 = points.c34(self.data.body, self._baseSystem, self.p33, self.p34)
        self.c35 = points.c35(self.data.body, self._baseSystem, self.p33, self.p34)
        self.c36, self.c37 = points.c36_37(self.p34, aux7, aux8, self.p36, t35)
        self.c38, self.c39 = points.c38_39(self.p34, aux7, aux8, self.p36, t35)
        self.c40 = points.c40(self.data.body, self._baseSystem, self.p36, self.p37)
        self.c41 = points.c41(self.data.body, self._baseSystem, self.p37, self.p37)

        self.c42e = points.c42(self.data.body, self._baseSystem, self.p1e, self.p29)
        self.c43e = points.c43(self.data.body, self._baseSystem, self.p1e, self.p29)
        self.c44e = points.c44(self.data.body, self._baseSystem, self.p1e, self.p34)
        self.c45e = points.c45(self.data.body, self._baseSystem, self.p1e, self.p34)
        self.c42d = self._symmetrical_point(self.c42e)
        self.c43d = self._symmetrical_point(self.c43e)
        self.c44d = self._symmetrical_point(self.c44e)
        self.c45d = self._symmetrical_point(self.c45e)

        self.c46e, self.c47e, self.c48e, self.c49e = points.c46_47_48_49(self.data.body, self._baseSystem, self.p30, self.p14e, self.p15e, self.p35, t30)
        self.c46d = self._symmetrical_point(self.c46e)
        self.c47d = self._symmetrical_point(self.c47e)
        self.c48d = self._symmetrical_point(self.c48e)
        self.c49d = self._symmetrical_point(self.c49e)

        self.c50e = points.c50(self.data.body, self._baseSystem, self.p13e, self.p31)
        self.c51e = points.c51(self.data.body, self._baseSystem, self.p13e, self.p31)
        self.c52e = points.c52(self.data.body, self._baseSystem, self.p13e, self.p36)
        self.c53e = points.c53(self.data.body, self._baseSystem, self.p13e, self.p36)
        self.c50d = self._symmetrical_point(self.c50e)
        self.c51d = self._symmetrical_point(self.c51e)
        self.c52d = self._symmetrical_point(self.c52e)
        self.c53d = self._symmetrical_point(self.c53e)
        
        return
    
    def _build_head(self) -> None:

        center = 0.5 * (self.p26e + self.p26d)
        self.p38 = points.p38(self.data.head, self._baseSystem, center)
        self.p39 = points.p39(self.data.head, self._baseSystem, center)
        self.p41e = points.p41e(self.data.head, self._baseSystem, center)
        self.p41d = points.p41d(self.data.head, self._baseSystem, center)
        self.p42 = points.p42(self.data.head, self._baseSystem, center)
        self.c54 = points.c54(self.data.head, self._baseSystem, self.p28)
        self.c55 = points.c55(self.data.head, self._baseSystem, self.p33)
        self.c56e = points.c56e(self.data.head, self._baseSystem, self.p26e)
        self.c56d = points.c56d(self.data.head, self._baseSystem, self.p26d)
        self.c57 = points.c57(self.data.head, self._baseSystem, self.p38)
        self.c58 = points.c58(self.data.head, self._baseSystem, self.p39)
        self.c59e = points.c59e(self.data.head, self._baseSystem, self.p41e)
        self.c59d = points.c59d(self.data.head, self._baseSystem, self.p41d)

        return
    
    def _build_tail(self) -> None:

        system = TailSystem(
            theta8=self.data.tail.theta8,
            theta9=self.data.tail.theta9,
            theta10=self.data.tail.theta10,
        )

        center = 0.5 * (self.p32 + self.p37)
        self.p47e = points.p47e(self.data.tail, system, center)
        self.p49e = points.p49e(self.data.body, self.data.tail, system, center)
        self.p51 = points.p51(self.data.tail, system, center)
        self.p49d = points.p49d(self.data.body, self.data.tail, system, center)
        self.p47d = points.p47d(self.data.tail, system, center)

        self.p43e, self.p44e, self.curve84e, self.curve85e, self.curve86e, self.curve87e = sections.process_section(self.data.tail.foil, self.p27e, self.p49e, -system.x4, system.z4, n=100)
        self.p43d, self.p44d, self.curve84d, self.curve85d, self.curve86d, self.curve87d = sections.process_section(self.data.tail.foil, self.p27d, self.p49d, -system.x4, system.z4, n=100)
        
        self.p45 = points.p45(self.p32, self.p51)
        self.p46 = points.p46(self.p37, self.p51)

        self.curve82e = points.curve82e(self.data.body, self.data.tail, system, center, 100)
        self.curve82d = points.curve82d(self.data.body, self.data.tail, system, center, 100)
        self.curve83e = points.curve83e(self.data.body, self.data.tail, system, center, 50)
        self.curve83d = points.curve83d(self.data.body, self.data.tail, system, center, 50)

        if self.data.tail.shape == TailShape.rounded:
            self.p51 = 0.5 * (self.curve83e[49, :] + self.curve83d[49, :])
        
        return