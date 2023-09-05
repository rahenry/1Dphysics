#!/usr/bin/env python3

class free_test_point:
    def __init__(self):
        self.sum = []
        1

def free_search(excitations, target, blocked=[]):
    tol = 1E-12
    if target.real < 0: return None
    for e, ind in excitations:
        if (abs(e) < tol):
            continue
        if ind in blocked:
            continue
        z = target-e
        if abs(z) < tol:
            return blocked+[ind]
        if (z.real < 0):
            break
        F = free_search(excitations, target-e, blocked+[ind])
        # F = free_search(excitations, target-e, [])
        if F:
            return F

    return None

def free_test_s(s, **kwargs):
    eigenstates = [s[f'e{i}'] for i in range(s['N'] ** s['L'])]
    return free_test_newnew(eigenstates, s['N'], **kwargs)

def free_test(z, N, use_excitations=False, base=0):
    print(f"Starting free_test, {base}")
    energies = sorted(z['energies'], key=lambda x:x.real)
    energies = energies[base:]
    # print(energies)

    e0 = energies[0]
    excitations = [e for e in energies[1:]]
    if use_excitations:
        excitations = [e - e0 for e in energies[1:]]
    # excitations = [e - e0 for e in energies[1:]]
    pfs = []
    pfs_dict = {}
    composites = {}
    ftps = []
    ind = 0
    res = []
    for e in excitations:
        print(e)
        sr = free_search(pfs, e, [])
        if sr and len(sr) > 1:
            s = sum([pfs_dict[i] for i in sr])
            # print(e, sr, s)
            res.append([e, sr])
        elif sr and len(sr) == 1:
            res.append([e, sr])
        else:
            # print('p', e)
            pfs.append([e, ind])
            pfs_dict[ind] = e
            if (abs(e) < 1E-6):
                res.append([e, []])
            else:
                res.append([e, [ind]])
            ind += 1
            # print([e, [ind]])

    def q(x):
        return [x[0], [y+1 for y in x[1]]]
    res = [q(x) for x in res]
    res = [[e0, [0]]] + res
    return res


def test1(e, others, stack=[], ignore=[]):
    tol = 1E-6
    tol2 = 1E-7
    q = e - sum(stack)
    for x in others:
        if x in ignore: continue
        # avoid dupes
        dupe = False
        if dupe: continue
        if (abs(q-x) < tol):
            return stack + [x]
        others_new = list(others)
        others_new.remove(x)
        others_new = [y for y in others_new if rel_diff(y, -x) > tol2]
        z = test1(e, others_new, stack+[x], list(ignore))
        if z:
            return z
        else:
            ignore.append(x)
    return None


def sum_test(z):
    N = z.N
    energies = sorted(z.eigenvalues, key=lambda x:x.real)
    e0 = energies[0]
    excitations = [e - e0 for e in energies[1:]]

    tol1 = 1E-8
    res = []
    energies = [x.real for x in energies if not abs(x) < tol1]
    for e in energies:
        others = [x for x in energies if not x==e]
        # dupes = []
        # for o in others:
        #     others2 = list(others)
        #     others2.remove(o)
        #     for v in others2:
        #         if rel_diff(o, -v) < tol1:
        #             dupes.append(o)
        #             dupes.append(v)
        # print(e, dupes)
        # continue
        qq = 1.000
        if rel_diff(e, qq) < 1E-8:
            t = test1(e, others)
            print(1)
            print(t)
        # continue
        t = test1(e, others)
        if t:
            res.append([e, t])
            print(t)
        # if t:
        #     print('...')
        #     print(e)
        #     print(t)
        #     print(sum(t))

    exit()
    report = ''
    report += f'L={z.L}, gamma={z.gamma}, bc={z.bc}\n'
    report += f'Eigenvalues:\n'
    for e in energies:
        report += f'{e.real:3.8f}\n'
    for (e, stack) in res:
        report += f'{e.real:3.8f} = {stack}\n'
    print(report)
    return res

def free_search_new(pfs, target, used=[], depth=0, start_ind=0):
    # print(depth, target)
    tol = 1E-8
    ind2 = 0
    for x, ind in pfs[start_ind:]:
        if abs(x) < 1E-6: continue
        if ind in used: continue
        d = target - x
        if abs(d) < tol:
            return used + [ind]
        if d.real < 0:
            break
        F = free_search_new(pfs, target-x, used+[ind], depth+1, ind2)
        if F:
            return F
    return None

def free_test_new(eigenstates, N, base=0):
    # energies = [e.energy for e in eigenstates]
    energies = eigenstates
    energies = sorted(energies, key=lambda x:x.real)
    energies = energies[base:]
    e0 = energies[0]
    excitations = [e - e0 for e in energies[1:]]

    ind = 0
    pfs_known = []
    pfs_dict = {}
    res = []
    for e in excitations:
        z = free_search_new(pfs_known, e, [])

        # This was a combination of other excitations (maybe just one)
        if z and len(z) >= 2:
            res.append([e, z])
        # This is a new excitation
        else:
            pf = [e, [ind]]
            pfs_known.append([e, ind])
            pfs_dict[ind] = e
            res.append(pf)
            ind += 1

    def q(x):
        return [x[0], [y+1 for y in x[1]]]
    res = [q(x) for x in res]
    # for i in range(len(res)):
    #     # s.measurement_data[f'free_data{i}'] = res[i]
    #     print(res[i])

    return res

def free_test_newnew(eigenstates, N, base=0):
    # energies = [e.energy for e in eigenstates]
    energies = eigenstates
    energies = sorted(energies, key=lambda x:x.real)
    energies = energies[base:]
    e0 = energies[0]
    excitations = [e - e0 for e in energies[1:]]

    tol = 1E-13
    pfs_known = []
    pfs_dict = {}
    res = []
    have_zero = False
    for e in excitations:
        if abs(e) < tol:
            have_zero = True

    zerod = {}
    by_e_ind = {}
    if have_zero:
        e_ind = 1
        while True:
            if e_ind >= len(excitations):
                break
            e = excitations[e_ind]
            e2_ind = e_ind+1
            found = False
            for e2 in excitations[e_ind+1:]:
                if (abs(e - e2) < tol):
                    zerod[e_ind] = e2_ind
                    e_ind = e2_ind+1
                    found = True
                    break

                e2_ind += 1

            if not found:
                e_ind += 1

    zerod_rev = {}
    for x, y in zerod.items():
        zerod_rev[y] = x
    e_ind = 0
    ind = 1
    for e in excitations:
        if e_ind in zerod_rev:
            print(f'Skipping {e_ind}')
            e_ind += 1
            continue
        if e_ind == 0 and have_zero:
            e_ind += 1
            continue

        z = free_search_new(pfs_known, e, [])

        # This was a combination of other excitations (maybe just one)
        if z and len(z) >= 2:
            res.append([e, z])
            by_e_ind[e_ind] = [e, z]
        # This is a new excitation
        else:
            pf = [e, [ind]]
            pfs_known.append([e, ind])
            pfs_dict[ind] = e
            res.append(pf)
            by_e_ind[e_ind] = pf
            ind += 1
        e_ind += 1

    if have_zero:
        res.append([0, [0]])
        for x, y in zerod.items():
            z = by_e_ind[x]
            res.append([z[0], [0]+z[1]])

    degens = {}
    d = 0
    seen = set()
    for p1 in pfs_known:
        i1 = p1[1]
        if i1 not in seen:
            seen.add(i1)
            degens[i1] = i1-d
        for p2 in pfs_known:
            i2 = p2[1]
            if not i2 >= i1: continue
            if p1[1] == p2[1]: continue
            if abs(p1[0] - p2[0]) < tol:
                degens[p2[1]] = p1[1] - d
                d += 1

    def replace_degens(x):
        res = []
        for y in x[1]:
            if y in degens:
                res.append(degens[y])
            else:
                res.append(y)
        return [x[0], res]
    res = [replace_degens(x) for x in res]
    # def q(x):
    #     return [x[0], [y+1 for y in x[1]]]
    # res = [q(x) for x in res]
    # for i in range(len(res)):
    #     # s.measurement_data[f'free_data{i}'] = res[i]
    #     print(res[i])

    return res

def free_test_reindex(eigs):
    1
