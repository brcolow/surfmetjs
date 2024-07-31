import numpy as np
import matplotlib.pyplot as plt

def cartesian_to_polar(points):
  """Converts an array of 2D Cartesian points to polar coordinates.

  Args:
    points: A NumPy array of shape (n, 2), where n is the number of points.

  Returns:
    A NumPy array of shape (n, 2), where the first column is the radius and the second column is the angle in radians.
  """

  x, y = points[:, 0], points[:, 1]
  r = np.sqrt(x**2 + y**2)
  theta = np.arctan2(y, x)
  return np.column_stack((r, theta))
    
n_samples = 500

nist_cir_data = np.array([[-37.8542, 1.94727, 2.6954],
[-37.58897, 1.88988, 2.6954],
[-37.40469, 1.86389, 2.6954],
[-37.54657, 1.48975, 2.6954],
[-37.19397, 1.5742, 2.6954],
[-37.00018, 1.52233, 2.6954],
[-37.07928, 1.15805, 2.6954],
[-36.78475, 1.25434, 2.6954],
[-36.82794, 0.97102, 2.6954],
[-36.48364, 0.98987, 2.6954],
[-36.24921, 0.98224, 2.6954],
[-36.38725, 0.54197, 2.6954],
[-36.23209, 0.46831, 2.6954],
[-35.94027, 0.5168, 2.6954],
[-35.70161, 0.51432, 2.6954],
[-35.80588, 0.09886, 2.6954],
[-35.48517, 0.23579, 2.6954],
[-35.36536, 0.12248, 2.6954],
[-35.26585, -0.0388, 2.6954],
[-34.96465, 0.09345, 2.6954],
[-34.89847, -0.18901, 2.6954],
[-34.86548, -0.53563, 2.6954],
[-34.68071, -0.42212, 2.6954],
[-34.44997, -0.54822, 2.6954],
[-34.32476, -0.58787, 2.6954],
[-34.0828, -0.63745, 2.6954],
[-33.89533, -0.64344, 2.6954],
[-33.90506, -0.98865, 2.6954],
[-33.90695, -1.30695, 2.6954],
[-33.73923, -1.40706, 2.6954],
[-33.32835, -1.12787, 2.6954],
[-33.1506, -1.1479, 2.6954],
[-33.01654, -1.22587, 2.6954],
[-32.95003, -1.48178, 2.6954],
[-32.78517, -1.62758, 2.6954],
[-32.52341, -1.49466, 2.6954],
[-32.53622, -1.96348, 2.6954],
[-32.28034, -1.8074, 2.6954],
[-32.29555, -2.20215, 2.6954],
[-32.01416, -2.06952, 2.6954],
[-31.82698, -2.16278, 2.6954],
[-31.73445, -2.3365, 2.6954],
[-31.53253, -2.31876, 2.6954],
[-31.32896, -2.33244, 2.6954],
[-31.12332, -2.37907, 2.6954],
[-31.0933, -2.72984, 2.6954],
[-30.79198, -2.50331, 2.6954],
[-30.70088, -2.65783, 2.6954],
[-30.65816, -2.95543, 2.6954],
[-30.33294, -2.70569, 2.6954],
[-30.15514, -2.80156, 2.6954],
[-29.98748, -2.88123, 2.6954],
[-29.9146, -3.17113, 2.6954],
[-29.71693, -3.16075, 2.6954],
[-29.58899, -3.23934, 2.6954],
[-29.35059, -3.2969, 2.6954],
[-29.22255, -3.4496, 2.6954],
[-29.16246, -3.60551, 2.6954],
[-28.84927, -3.4423, 2.6954],
[-28.74377, -3.73398, 2.6954],
[-28.51492, -3.70852, 2.6954],
[-28.33995, -3.6943, 2.6954],
[-28.18017, -3.68529, 2.6954],
[-27.96923, -3.72787, 2.6954],
[-27.92217, -4.07808, 2.6954],
[-27.62104, -3.77041, 2.6954],
[-27.62531, -4.26556, 2.6954],
[-27.25452, -3.83936, 2.6954],
[-27.20122, -4.27652, 2.6954],
[-27.0374, -4.22096, 2.6954],
[-26.72695, -3.9724, 2.6954],
[-26.57394, -4.02114, 2.6954],
[-26.50773, -4.48666, 2.6954],
[-26.19798, -4.21668, 2.6954],
[-26.12512, -4.4042, 2.6954],
[-25.89375, -4.30562, 2.6954],
[-25.7199, -4.32215, 2.6954],
[-25.64268, -4.75042, 2.6954],
[-25.37149, -4.43433, 2.6954],
[-25.21085, -4.35692, 2.6954],
[-25.06462, -4.72192, 2.6954],
[-24.88036, -4.50221, 2.6954],
[-24.67389, -4.75544, 2.6954],
[-24.57363, -4.96182, 2.6954],
[-24.3835, -5.03767, 2.6954],
[-24.14689, -4.78337, 2.6954],
[-24.0146, -4.99425, 2.6954],
[-23.8514, -5.18314, 2.6954],
[-23.56931, -4.71068, 2.6954],
[-23.4925, -5.02922, 2.6954],
[-23.29964, -4.96341, 2.6954],
[-23.07176, -5.01766, 2.6954],
[-22.8913, -5.29753, 2.6954],
[-22.7218, -5.00778, 2.6954],
[-22.54225, -4.98807, 2.6954],
[-22.38065, -5.16108, 2.6954],
[-22.14858, -5.0746, 2.6954],
[-22.0507, -5.35625, 2.6954],
[-21.84989, -5.46464, 2.6954],
[-21.69893, -5.48884, 2.6954],
[-21.50466, -5.22397, 2.6954],
[-21.25848, -5.19547, 2.6954],
[-21.10176, -5.11473, 2.6954],
[-20.95506, -5.52511, 2.6954],
[-20.71586, -5.24515, 2.6954],
[-20.58004, -5.34917, 2.6954],
[-20.3809, -5.14646, 2.6954],
[-20.21413, -5.18603, 2.6954],
[-19.95904, -5.1767, 2.6954],
[-19.79424, -5.65261, 2.6954],
[-19.66336, -5.10619, 2.6954],
[-19.48105, -5.61969, 2.6954],
[-19.27754, -5.19302, 2.6954],
[-19.10153, -5.48969, 2.6954],
[-18.91209, -5.64348, 2.6954],
[-18.74852, -5.2085, 2.6954],
[-18.56682, -5.55356, 2.6954],
[-18.40054, -5.30175, 2.6954],
[-18.20079, -5.58949, 2.6954],
[-17.99361, -5.60221, 2.6954],
[-17.84044, -5.34374, 2.6954],
[-17.65064, -5.56322, 2.6954],
[-17.49352, -5.13811, 2.6954],
[-17.25381, -5.54983, 2.6954],
[-17.11511, -5.45764, 2.6954],
[-16.93388, -5.25844, 2.6954],
[-16.77096, -5.1302, 2.6954],
[-16.56586, -5.19889, 2.6954],
[-16.36595, -5.61585, 2.6954],
[-16.20395, -5.11801, 2.6954],
[-15.99513, -5.50775, 2.6954],
[-15.86288, -5.56603, 2.6954],
[-15.66119, -5.18766, 2.6954],
[-15.45772, -5.07597, 2.6954],
[-15.34884, -5.22085, 2.6954],
[-15.11098, -4.99589, 2.6954],
[-14.9869, -5.06531, 2.6954],
[-14.74444, -5.32393, 2.6954],
[-14.53963, -5.46995, 2.6954],
[-14.32389, -5.39686, 2.6954],
[-14.20332, -5.16567, 2.6954],
[-13.97246, -5.31445, 2.6954],
[-13.84589, -5.34296, 2.6954],
[-13.64285, -5.19185, 2.6954],
[-13.55387, -4.73014, 2.6954],
[-13.33269, -4.91355, 2.6954],
[-13.11383, -5.08386, 2.6954],
[-12.96737, -4.87898, 2.6954],
[-12.74999, -5.13455, 2.6954],
[-12.62602, -4.73729, 2.6954],
[-12.50452, -4.67084, 2.6954],
[-12.27001, -4.63853, 2.6954],
[-12.11556, -4.42056, 2.6954],
[-11.97239, -4.54118, 2.6954],
[-11.72314, -4.55026, 2.6954],
[-11.55783, -4.58693, 2.6954],
[-11.33432, -4.62624, 2.6954],
[-11.13842, -4.6212, 2.6954],
[-10.94379, -4.7409, 2.6954],
[-10.87307, -4.13824, 2.6954],
[-10.68918, -4.3733, 2.6954],
[-10.46333, -4.52587, 2.6954],
[-10.36146, -4.096, 2.6954],
[-10.12685, -4.13876, 2.6954],
[-9.87363, -4.36652, 2.6954],
[-9.80715, -4.21472, 2.6954],
[-9.61908, -4.13735, 2.6954],
[-9.37188, -4.16839, 2.6954],
[-9.26289, -3.97754, 2.6954],
[-9.06989, -3.95748, 2.6954],
[-8.97204, -3.66792, 2.6954],
[-8.66657, -3.92863, 2.6954],
[-8.63993, -3.46867, 2.6954],
[-8.4676, -3.51935, 2.6954],
[-8.17875, -3.69254, 2.6954],
[-8.10416, -3.51723, 2.6954],
[-7.97189, -3.2768, 2.6954],
[-7.74795, -3.26614, 2.6954],
[-7.56937, -3.30887, 2.6954],
[-7.34412, -3.3209, 2.6954],
[-7.18087, -3.32471, 2.6954],
[-7.12666, -2.92343, 2.6954],
[-6.8805, -3.07555, 2.6954],
[-6.72534, -2.96729, 2.6954],
[-6.47907, -3.17591, 2.6954],
[-6.4943, -2.65754, 2.6954],
[-6.17563, -2.90865, 2.6954],
[-6.11296, -2.601, 2.6954],
[-5.95338, -2.55642, 2.6954],
[-5.75482, -2.48342, 2.6954],
[-5.62691, -2.38691, 2.6954],
[-5.42179, -2.31226, 2.6954],
[-5.17782, -2.41502, 2.6954],
[-5.09443, -2.10808, 2.6954],
[-4.93696, -2.07278, 2.6954],
[-4.90484, -1.76174, 2.6954],
[-4.62638, -1.77313, 2.6954],
[-4.54943, -1.64622, 2.6954],
[-4.38112, -1.62499, 2.6954],
[-4.28002, -1.39576, 2.6954],
[-4.06802, -1.37891, 2.6954],
[-3.85284, -1.36721, 2.6954],
[-3.54542, -1.57082, 2.6954],
[-3.4066, -1.4057, 2.6954],
[-3.32674, -1.27862, 2.6954],
[-3.23725, -1.03399, 2.6954],
[-3.14691, -0.85603, 2.6954],
[-2.93798, -0.83438, 2.6954],
[-2.81187, -0.71604, 2.6954],
[-2.72932, -0.44653, 2.6954],
[-2.51581, -0.53537, 2.6954],
[-2.30685, -0.49969, 2.6954],
[-2.23796, -0.18751, 2.6954],
[-1.9574, -0.38117, 2.6954],
[-1.88814, -0.11743, 2.6954],
[-1.56595, -0.301, 2.6954],
[-1.62132, 0.14669, 2.6954],
[-1.54351, 0.34883, 2.6954],
[-1.32749, 0.32944, 2.6954],
[-0.98308, 0.25832, 2.6954],
[-1.12222, 0.71119, 2.6954],
[-0.89244, 0.69835, 2.6954],
[-0.64143, 0.69512, 2.6954],
[-0.60178, 0.90307, 2.6954],
[-0.47029, 1.01585, 2.6954],
[-0.43434, 1.2349, 2.6954],
[-0.17733, 1.23361, 2.6954],
[-0.07842, 1.39041, 2.6954],
[0.30197, 1.27154, 2.6954],
[0.33127, 1.50287, 2.6954],
[0.39671, 1.67559, 2.6954],
[0.72226, 1.61975, 2.6954],
[0.85757, 1.78344, 2.6954],
[0.79516, 2.04351, 2.6954],
[0.96201, 2.12186, 2.6954],
[1.30235, 2.08241, 2.6954],
[1.37888, 2.16561, 2.6954],
[1.35561, 2.53848, 2.6954],
[1.61696, 2.47622, 2.6954],
[1.7001, 2.6712, 2.6954],
[1.60515, 3.03066, 2.6954],
[1.76284, 3.15023, 2.6954],
[1.9956, 3.08717, 2.6954],
[2.01852, 3.37002, 2.6954],
[2.45635, 3.22708, 2.6954],
[2.35827, 3.58909, 2.6954],
[2.38984, 3.70489, 2.6954],
[2.83319, 3.60166, 2.6954],
[2.96478, 3.73831, 2.6954],
[2.76767, 4.16681, 2.6954],
[2.99393, 4.20861, 2.6954],
[3.11186, 4.34645, 2.6954],
[3.40229, 4.32457, 2.6954],
[3.45366, 4.46839, 2.6954],
[3.45743, 4.76888, 2.6954],
[3.7735, 4.7432, 2.6954],
[3.53155, 5.14688, 2.6954],
[4.07694, 5.00585, 2.6954],
[3.82086, 5.38438, 2.6954],
[4.09235, 5.35699, 2.6954],
[4.13137, 5.59187, 2.6954],
[4.21882, 5.76303, 2.6954],
[4.58515, 5.69234, 2.6954],
[4.65298, 5.91042, 2.6954],
[4.76117, 6.06914, 2.6954],
[4.54332, 6.40478, 2.6954],
[4.92909, 6.45248, 2.6954],
[4.78812, 6.72065, 2.6954],
[5.05026, 6.74235, 2.6954],
[5.24226, 6.82952, 2.6954],
[5.39853, 6.91266, 2.6954],
[5.14028, 7.35794, 2.6954],
[5.62079, 7.17712, 2.6954],
[5.43103, 7.58962, 2.6954],
[5.68189, 7.57225, 2.6954],
[5.94619, 7.68589, 2.6954],
[5.73818, 8.07012, 2.6954],
[5.91363, 8.12776, 2.6954],
[5.97114, 8.3207, 2.6954],
[6.3687, 8.30597, 2.6954],
[6.29522, 8.52781, 2.6954],
[6.4003, 8.69335, 2.6954],
[6.3773, 8.87799, 2.6954],
[6.44847, 9.03591, 2.6954],
[6.47944, 9.2657, 2.6954],
[6.73967, 9.33934, 2.6954],
[6.78685, 9.45026, 2.6954],
[6.88577, 9.65756, 2.6954],
[6.77129, 9.92423, 2.6954],
[7.23786, 9.85305, 2.6954],
[7.34746, 10.0341, 2.6954],
[7.12593, 10.34639, 2.6954],
[7.12398, 10.57348, 2.6954],
[7.59323, 10.52929, 2.6954],
[7.57768, 10.72061, 2.6954],
[7.70877, 10.86544, 2.6954],
[7.49463, 11.21668, 2.6954],
[7.55917, 11.35553, 2.6954],
[7.99836, 11.35304, 2.6954],
[7.71754, 11.66444, 2.6954],
[7.83727, 11.81768, 2.6954],
[7.92601, 11.98775, 2.6954],
[7.98696, 12.15084, 2.6954],
[8.31718, 12.18408, 2.6954],
[8.20835, 12.39948, 2.6954],
[8.2597, 12.60107, 2.6954],
[8.50184, 12.72947, 2.6954],
[8.51803, 12.91868, 2.6954],
[8.67676, 13.07052, 2.6954],
[8.32359, 13.39747, 2.6954],
[8.63593, 13.43729, 2.6954],
[8.66148, 13.58518, 2.6954],
[8.77075, 13.79291, 2.6954],
[8.96905, 13.8895, 2.6954],
[9.00513, 14.1114, 2.6954],
[8.8294, 14.32799, 2.6954],
[8.8546, 14.55246, 2.6954],
[9.08692, 14.65969, 2.6954],
[8.83112, 14.9578, 2.6954],
[9.39485, 14.91184, 2.6954],
[9.10591, 15.26289, 2.6954],
[9.05297, 15.43064, 2.6954],
[9.0685, 15.59023, 2.6954],
[9.59052, 15.67291, 2.6954],
[9.61126, 15.85105, 2.6954],
[9.26543, 16.11163, 2.6954],
[9.73138, 16.1651, 2.6954],
[9.83719, 16.35614, 2.6954],
[9.76482, 16.59488, 2.6954],
[9.45081, 16.80142, 2.6954],
[9.81785, 16.88121, 2.6954],
[9.8822, 17.09814, 2.6954],
[9.9548, 17.2616, 2.6954],
[10.04501, 17.41974, 2.6954],
[10.04635, 17.59631, 2.6954],
[9.74446, 17.85783, 2.6954],
[9.82256, 18.03337, 2.6954],
[9.83072, 18.21156, 2.6954],
[10.06879, 18.37529, 2.6954],
[9.84359, 18.59208, 2.6954],
[10.02441, 18.7685, 2.6954],
[10.27173, 18.85736, 2.6954],
[9.90457, 19.10873, 2.6954],
[9.89563, 19.32175, 2.6954],
[10.04511, 19.50085, 2.6954],
[10.28407, 19.65071, 2.6954],
[10.24855, 19.82222, 2.6954],
[10.01736, 20.00628, 2.6954],
[10.30878, 20.15822, 2.6954],
[10.179, 20.31056, 2.6954],
[10.07714, 20.54614, 2.6954],
[10.2913, 20.74595, 2.6954],
[10.29389, 20.89513, 2.6954],
[10.50301, 21.03082, 2.6954],
[10.42838, 21.26197, 2.6954],
[10.25429, 21.41946, 2.6954],
[10.57934, 21.63235, 2.6954],
[10.07864, 21.80895, 2.6954],
[10.21196, 21.96058, 2.6954],
[10.31773, 22.13933, 2.6954],
[10.27971, 22.34172, 2.6954],
[10.63653, 22.47639, 2.6954],
[10.12223, 22.6971, 2.6954],
[10.36614, 22.89953, 2.6954],
[10.42741, 23.08153, 2.6954],
[10.19985, 23.22683, 2.6954],
[10.61357, 23.43469, 2.6954],
[10.50699, 23.57301, 2.6954],
[10.25028, 23.81011, 2.6954],
[10.42108, 24.01289, 2.6954],
[10.46949, 24.19706, 2.6954],
[10.63584, 24.34619, 2.6954],
[10.43548, 24.5751, 2.6954],
[10.3185, 24.69607, 2.6954],
[10.39633, 24.89288, 2.6954],
[10.44597, 25.05015, 2.6954],
[10.24516, 25.2295, 2.6954],
[10.33431, 25.43477, 2.6954],
[10.42268, 25.63963, 2.6954],
[10.403, 25.81102, 2.6954],
[10.39172, 25.96036, 2.6954],
[10.56319, 26.21587, 2.6954],
[9.96747, 26.29702, 2.6954],
[10.10922, 26.48312, 2.6954],
[10.28043, 26.70934, 2.6954],
[10.01246, 26.83618, 2.6954],
[10.22057, 27.02802, 2.6954],
[9.87846, 27.22605, 2.6954],
[10.10007, 27.42894, 2.6954],
[9.97082, 27.60521, 2.6954],
[10.02185, 27.78041, 2.6954],
[10.27611, 28.03579, 2.6954],
[10.16572, 28.16306, 2.6954],
[10.02821, 28.32161, 2.6954],
[10.23758, 28.54468, 2.6954],
[10.09376, 28.74017, 2.6954],
[9.84914, 28.82303, 2.6954],
[9.80151, 29.0155, 2.6954],
[9.71404, 29.17376, 2.6954],
[9.58705, 29.29373, 2.6954],
[9.74396, 29.57387, 2.6954],
[9.78724, 29.79564, 2.6954],
[9.43758, 29.82291, 2.6954],
[9.61783, 30.08814, 2.6954],
[9.49842, 30.24345, 2.6954],
[9.33962, 30.43388, 2.6954],
[9.77467, 30.70943, 2.6954],
[9.27707, 30.7644, 2.6954],
[9.50111, 30.98178, 2.6954],
[9.30563, 31.11515, 2.6954],
[9.13336, 31.26478, 2.6954],
[9.29052, 31.48475, 2.6954],
[9.12335, 31.5957, 2.6954],
[8.93094, 31.78347, 2.6954],
[9.23123, 32.02066, 2.6954],
[8.84841, 32.12062, 2.6954],
[9.1253, 32.41783, 2.6954],
[8.77154, 32.50293, 2.6954],
[8.97197, 32.74367, 2.6954],
[9.01159, 32.98041, 2.6954],
[8.47888, 32.94924, 2.6954],
[8.96638, 33.33022, 2.6954],
[8.31888, 33.29604, 2.6954],
[8.4784, 33.48594, 2.6954],
[8.52023, 33.76793, 2.6954],
[8.50483, 33.93692, 2.6954],
[8.20729, 34.01218, 2.6954],
[8.23931, 34.26307, 2.6954],
[8.12397, 34.36815, 2.6954],
[8.29021, 34.60409, 2.6954],
[8.09545, 34.74788, 2.6954],
[8.00378, 34.92554, 2.6954],
[8.03182, 35.1154, 2.6954],
[8.0019, 35.3078, 2.6954],
[7.82942, 35.47416, 2.6954],
[7.64328, 35.57276, 2.6954],
[7.85609, 35.86665, 2.6954],
[7.70124, 36.00883, 2.6954],
[7.27073, 35.95012, 2.6954],
[7.40446, 36.23314, 2.6954],
[7.37414, 36.392, 2.6954],
[7.23125, 36.56025, 2.6954],
[7.39104, 36.79948, 2.6954],
[7.03393, 36.85127, 2.6954],
[7.25134, 37.19354, 2.6954],
[7.08762, 37.27468, 2.6954],
[7.01693, 37.41981, 2.6954],
[6.72363, 37.53348, 2.6954],
[6.3669, 37.55896, 2.6954],
[6.40372, 37.71475, 2.6954],
[6.43485, 37.99445, 2.6954],
[6.55659, 38.22973, 2.6954],
[6.02005, 38.13854, 2.6954],
[6.12184, 38.44553, 2.6954],
[6.02006, 38.6028, 2.6954],
[6.1646, 38.82591, 2.6954],
[5.86475, 38.93753, 2.6954],
[5.63129, 38.93694, 2.6954],
[5.814, 39.35437, 2.6954],
[5.80132, 39.54123, 2.6954],
[5.61574, 39.59116, 2.6954],
[5.58548, 39.81653, 2.6954],
[5.225, 39.77198, 2.6954],
[5.16537, 39.95106, 2.6954],
[4.89727, 39.99975, 2.6954],
[5.15413, 40.41468, 2.6954],
[4.73409, 40.33237, 2.6954],
[4.79412, 40.59828, 2.6954],
[4.54825, 40.65656, 2.6954],
[4.72847, 40.98632, 2.6954],
[4.46188, 40.99014, 2.6954],
[4.34568, 41.19984, 2.6954],
[4.21221, 41.34577, 2.6954],
[4.15364, 41.48377, 2.6954],
[4.0109, 41.59267, 2.6954],
[4.0611, 41.87798, 2.6954],
[3.66166, 41.74063, 2.6954],
[3.44491, 41.87517, 2.6954],
[3.68852, 42.24656, 2.6954],
[3.59352, 42.43641, 2.6954],
[3.53186, 42.60404, 2.6954],
[3.1505, 42.5446, 2.6954],
[2.87205, 42.51878, 2.6954],
[2.73775, 42.6603, 2.6954],
[2.59445, 42.78998, 2.6954],
[2.56275, 43.04578, 2.6954],
[2.60386, 43.22994, 2.6954],
[2.5379, 43.52104, 2.6954],
[2.27372, 43.47783, 2.6954],
[2.09954, 43.54439, 2.6954],
[2.21888, 43.96038, 2.6954],
[2.08726, 44.05534, 2.6954],
[1.77119, 44.03189, 2.6954],
[1.59574, 44.12822, 2.6954],
[1.69597, 44.48955, 2.6954],
[1.5621, 44.55156, 2.6954],
[1.23826, 44.51048, 2.6954],
[0.96843, 44.42778, 2.6954],
[0.97381, 44.78267, 2.6954],
[1.08561, 45.11837, 2.6954]])

xy_nist_data = nist_cir_data[:, :2]
polar_nist_data = cartesian_to_polar(xy_nist_data)
radius_nist_data = polar_nist_data[:, :1]
print(radius_nist_data)
np_fft = np.fft.rfft2(radius_nist_data)
amplitudes = (2 / n_samples) * np.abs(np_fft) 
frequencies = np.fft.fftfreq(n_samples) * n_samples * 1
plt.subplots(subplot_kw={'projection': 'polar'})
#plt.plot(polar_nist_data)
plt.plot(frequencies[:len(frequencies) // 2], amplitudes[:len(np_fft) // 2])
plt.show()
plt.savefig("output_plot.png")
print(amplitudes)
