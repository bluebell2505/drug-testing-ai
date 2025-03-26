
import React from 'react';
import Navbar from '@/components/Navbar';
import Sidebar from '@/components/Sidebar';
import Footer from '@/components/Footer';
import { BarChart2, TrendingUp, Atom, Activity } from 'lucide-react';

const Analytics: React.FC = () => {
  return (
    <div className="flex flex-col min-h-screen bg-gray-900">
      <Navbar />
      
      <div className="flex flex-1">
        <Sidebar />
        
        <main className="flex-1 p-6">
          <h1 className="text-white text-3xl font-bold mb-6 flex items-center gap-2">
            <BarChart2 className="text-blue-400" size={28} />
            Analytics
          </h1>
          
          <div className="grid grid-cols-1 md:grid-cols-2 gap-6 mb-6">
            <div className="bg-gray-800/50 border border-gray-700 p-6 rounded-xl">
              <h2 className="text-white text-xl font-bold mb-4 flex items-center gap-2">
                <TrendingUp size={20} className="text-blue-400" />
                AI Model Performance
              </h2>
              <p className="text-gray-300">Detailed insights on molecule generation, model accuracy, and prediction trends coming soon.</p>
              
              <div className="mt-4 grid grid-cols-2 gap-4">
                <div className="bg-gray-700/30 p-3 rounded-lg">
                  <div className="text-gray-400 text-sm">Accuracy</div>
                  <div className="text-blue-400 text-xl font-bold">94.2%</div>
                </div>
                <div className="bg-gray-700/30 p-3 rounded-lg">
                  <div className="text-gray-400 text-sm">Response Time</div>
                  <div className="text-blue-400 text-xl font-bold">1.2s</div>
                </div>
              </div>
            </div>
            
            <div className="bg-gray-800/50 border border-gray-700 p-6 rounded-xl">
              <h2 className="text-white text-xl font-bold mb-4 flex items-center gap-2">
                <Atom size={20} className="text-blue-400" />
                Molecule Generation Stats
              </h2>
              <p className="text-gray-300">Statistics on generated molecules and their properties.</p>
              
              <div className="mt-4 grid grid-cols-2 gap-4">
                <div className="bg-gray-700/30 p-3 rounded-lg">
                  <div className="text-gray-400 text-sm">Total Generated</div>
                  <div className="text-blue-400 text-xl font-bold">128</div>
                </div>
                <div className="bg-gray-700/30 p-3 rounded-lg">
                  <div className="text-gray-400 text-sm">Success Rate</div>
                  <div className="text-blue-400 text-xl font-bold">87.5%</div>
                </div>
              </div>
            </div>
          </div>
          
          <div className="bg-gray-800/50 border border-gray-700 p-6 rounded-xl">
            <h2 className="text-white text-xl font-bold mb-4 flex items-center gap-2">
              <Activity size={20} className="text-blue-400" />
              Activity Timeline
            </h2>
            <p className="text-gray-300">
              The analytics dashboard is under development. Charts and detailed metrics on molecule generation, 
              model performance, and prediction trends will be available soon.
            </p>
          </div>
        </main>
      </div>
      
      <Footer />
    </div>
  );
};

export default Analytics;
