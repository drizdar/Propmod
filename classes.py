import formulas as f
class solution:
    def __init__(self, data):
        self.data = data
    def CalcFlow(self):
        data = self.data
        rho = data.get("density")
        m_t = data.get("mass_total")
        if (rho != None and m_t != None):
            data["flow"] = m_t / rho
            return True
        else:
            print(f'density = {rho}, m_t = {m_t}: could not complete calculation')
            return None
    def CalcPcWt(self):
        data = self.data
        pc_wt = data.get("pc_wt")
        if pc_wt != None:
            print("pc_wt already exists")
            return None
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
            return True
        else:
            print(f'mass_water = {m_w}, mass_NaCl = {m_NaCl}: could not complete calculation')
            return None
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
            return True
        else:
            print(f'flow = {flow}, density = {rho}, pc_wt = {pc_wt}: could not compelte calculation')
            return None

    def CalcOsmoticProperties(self):
        data = self.data
        P = data.get("P")
        T = data.get("T")
        pc_wt = data.get("pc_wt")
        if (P != None and T != None and pc_wt != None):
            data["PI"], data["density"] = f.OsmoticProperties(P, T, pc_wt)
            return True
        else:
            print(f'P = {P}, T = {T}, pc_wt = {pc_wt}: could not compelte calculation')
            return None