import formulas as f
import math
import Reporting as r
class solution:
    def __init__(self, data):
        self.data = data
    def CalcFlow(self):
        data = self.data
        rho = data.get("density")
        m_t = data.get("mass_total")
        if (rho != None and m_t != None):
            data["flow"] = m_t / rho
        else:
            raise Exception(f'density = {rho}, m_t = {m_t}: could not complete calculation')
    def CalcPcWt(self):
        data = self.data
        m_w = data.get("mass_water")
        m_NaCl = data.get("mass_NaCl")
        if (m_w != None and m_NaCl != None):
            m_t = m_w + m_NaCl
            data["mass_total"] = m_t
            try:                
                pc_wt = m_NaCl/m_t
            except:
                pc_wt = 0
            data["pc_wt"] = pc_wt
        else:
            raise Exception(f'mass_water = {m_w}, mass_NaCl = {m_NaCl}: could not complete calculation')
    def CalcMass(self):
        data = self.data
        pc_wt = data.get("pc_wt")
        rho = data.get("density")
        flow = data.get("flow")
        if (pc_wt != None and rho != None and flow != None):
            m_w = flow * rho * (1-pc_wt) * 1e6
            m_NaCl =  flow * rho * (pc_wt) * 1e6
            m_t = m_w + m_NaCl
            data["mass_water"] = m_w
            data["mass_NaCl"] = m_NaCl
            data["mass_total"] = m_t
        else:
            raise Exception(f'flow = {flow}, density = {rho}, pc_wt = {pc_wt}: could not complete calculation')

    def CalcOsmoticProperties(self):
        data = self.data
        P = data.get("P")
        T = data.get("T")
        pc_wt = data.get("pc_wt")
        if (P != None and T != None and pc_wt != None):
            data["PI"], data["density"] = f.OsmoticProperties(P, T, pc_wt)
        else:
            raise Exception(f'P = {P}, T = {T}, pc_wt = {pc_wt}: could not complete calculation')

class membrane:
    def __init__(self, data, s):
        self.data = data
        self.CalcLf()
        self.data["dimensions"]["s"] = s
        self.CalcNed()
        self.CalcNef()
        self.CalcDLd()
        self.CalcDLf()
        self.CalcDAd()
        self.CalcDAf()
        self.CalcDAm()
    def CalcD_h(self):
        dims = self.data["dimensions"]
        Hc = dims.get("Hc")
        Ss = dims.get("Ss")
        if (Hc != None and Ss != None):
            dims["d_h"] = f.HydraulicDiameter(Hc, Ss)
        else:
            raise Exception(f'Hc = {Hc}, SS = {Ss}: could not complete calculation')
    def CalcLf(self):
        dims = self.data["dimensions"]
        Am = dims.get("Am")
        Ld = dims.get("Ld")
        if (Am != None and Ld != None):
            dims["Lf"] = Am/Ld
        else:
            raise Exception(f'Am = {Am}, Ld = {Ld}: could not complete calculations')
    def CalcNed(self):
        dims = self.data["dimensions"]
        Ld = dims.get("Ld")
        s = dims.get("s")
        if (Ld != None and s != None):
            dims["ned"] = math.ceil(Ld/s)
        else:
            raise Exception(f'Ld = {Ld}, s = {s}: could not complete calculations')
    def CalcNef(self):
        dims = self.data["dimensions"]
        Ld = dims.get("Lf")
        s = dims.get("s")
        if (Ld != None and s != None):
            dims["nef"] = math.ceil(Ld/s)
        else:
            raise Exception(f'Ld = {Ld}, s = {s}: could not complete calculations')
    def CalcDLd(self):
        dims = self.data["dimensions"]
        Ld = dims.get("Ld")
        ned = dims.get("ned")
        if (Ld != None and ned != None):
            dims["dLd"] = Ld/ned
        else:
            raise Exception(f'Ld = {Ld}, ned = {ned}: could not complete calculations')
    def CalcDLf(self):
        dims = self.data["dimensions"]
        Lf = dims.get("Lf")
        nef = dims.get("nef")
        if (Lf != None and nef != None):
            dims["dLf"] = Lf/nef
        else:
            raise Exception(f'Ld = {Lf}, nef = {nef}: could not complete calculations')
    def CalcDAd(self):
        dims = self.data["dimensions"]
        Hc = dims.get("Hc")
        dLd = dims.get("dLd")
        if (Hc != None and dLd != None):
            dims["dAd"] = Hc*dLd
        else:
            raise Exception(f'Hc = {Hc}, dLd = {dLd}: could not complete calculations')
    def CalcDAf(self):
        dims = self.data["dimensions"]
        Hc = dims.get("Hc")
        dLf = dims.get("dLf")
        if (Hc != None and dLf != None):
            dims["dAf"] = Hc*dLf
        else:
            raise Exception(f'Hc = {Hc}, dLd = {dLf}: could not complete calculations')
    def CalcDAm(self):
        dims = self.data["dimensions"]
        dLd = dims.get("dLd")
        dLf = dims.get("dLf")
        if (dLd != None and dLf != None):
            dims["dAm"] = dLd*dLf
        else:
            raise Exception(f'dLd = {dLd}, dLd = {dLf}: could not complete calculations')
    def GetProperties(self, list):
        ret_array = []
        for item in list:
            if item in self.data["properties"]:
                ret_array.append(self.data["properties"][item])
        return ret_array
    def GetDimensions(self, list):
        ret_array = []
        for item in list:
            if (item in self.data["dimensions"] and type(self.data["dimensions"][item]) == "float"):
                ret_array.append(self.data["dimensions"][item])
            else:
                raise Exception(f'Item {item} was either missing or not of type "float"')
        return ret_array
    def Report(self, rep):
        data = self.data
        rep = r.pushReport(data["properties"]["A"], "A: Water permeability coefficient (L/h/m^2/bar)", rep)
        rep = r.pushReport(data["properties"]["B"], "B: Salt permeability coefficient (g/h/m^2)", rep)
        rep = r.pushReport(data["properties"]["S"], "S: Structural parameter", rep)
        rep = r.pushReport(data["dimensions"]["Am"], "Area of membrane (m^2)", rep)
        rep = r.pushReport(data["dimensions"]["Ld"], "Length of membrane along draw direction (m)", rep)
        rep = r.pushReport(data["dimensions"]["Lf"], "Length of membrane along feed direction (m)", rep)
        rep = r.pushReport(data["dimensions"]["Hc"],'H: Channel height (m)', rep)
        rep = r.pushReport(data["dimensions"]["Ss"],'Ss: Spacer Spacing (m)', rep)
        rep = r.pushReport(data["dimensions"]["s"], "Approximate length of calculation element (m)", rep)
        rep = r.pushReport(data["dimensions"]["ned"], "Number of elements in draw direction", rep)
        rep = r.pushReport(data["dimensions"]["nef"], "Number of elements in feed direction", rep)
        rep = r.pushReport(data["dimensions"]["dLd"], "Length of element in draw direction (m)", rep)
        rep = r.pushReport(data["dimensions"]["dLf"], "Length of element in feed direction (m)", rep)
        rep = r.pushReport(data["dimensions"]["dAd"], "dAd: Area of channel in draw direction (m^2)", rep)
        rep = r.pushReport(data["dimensions"]["dAf"], "dAf: Area of channel in feed direction (m^2)", rep)
        rep = r.pushReport(data["dimensions"]["dAm"], "dAm: Area of element (m^2)", rep)
        return rep