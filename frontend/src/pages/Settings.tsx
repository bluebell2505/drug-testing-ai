
import React from 'react';
import Navbar from '@/components/Navbar';
import Sidebar from '@/components/Sidebar';
import Footer from '@/components/Footer';
import { Settings as SettingsIcon, User, Database, Bell, Shield, Info } from 'lucide-react';
import { Button } from '@/components/ui/button';

const Settings: React.FC = () => {
  return (
    <div className="flex flex-col min-h-screen bg-gray-900">
      <Navbar />
      
      <div className="flex flex-1">
        <Sidebar />
        
        <main className="flex-1 p-6">
          <h1 className="text-white text-3xl font-bold mb-6 flex items-center gap-2">
            <SettingsIcon className="text-blue-400" size={28} />
            Settings
          </h1>
          
          <div className="grid grid-cols-1 lg:grid-cols-3 gap-6">
            <div className="lg:col-span-1">
              <div className="bg-gray-800/50 border border-gray-700 p-4 rounded-xl sticky top-6">
                <nav className="flex flex-col space-y-1">
                  <a href="#user" className="flex items-center gap-2 text-blue-400 bg-gray-700/30 p-3 rounded-lg font-medium">
                    <User size={18} />
                    User Settings
                  </a>
                  <a href="#api" className="flex items-center gap-2 text-gray-300 hover:text-blue-400 p-3 rounded-lg">
                    <Database size={18} />
                    API Configuration
                  </a>
                  <a href="#notifications" className="flex items-center gap-2 text-gray-300 hover:text-blue-400 p-3 rounded-lg">
                    <Bell size={18} />
                    Notifications
                  </a>
                  <a href="#privacy" className="flex items-center gap-2 text-gray-300 hover:text-blue-400 p-3 rounded-lg">
                    <Shield size={18} />
                    Privacy & Security
                  </a>
                  <a href="#about" className="flex items-center gap-2 text-gray-300 hover:text-blue-400 p-3 rounded-lg">
                    <Info size={18} />
                    About
                  </a>
                </nav>
              </div>
            </div>
            
            <div className="lg:col-span-2">
              <div id="user" className="bg-gray-800/50 border border-gray-700 p-6 rounded-xl mb-6">
                <h2 className="text-white text-2xl font-bold mb-4">User Settings</h2>
                
                <div className="space-y-4">
                  <div>
                    <label className="block text-white font-medium mb-2">Full Name</label>
                    <input 
                      type="text" 
                      value="Dr. John Smith"
                      className="w-full bg-gray-700 border border-gray-600 text-white p-3 rounded-lg"
                    />
                  </div>
                  
                  <div>
                    <label className="block text-white font-medium mb-2">Email</label>
                    <input 
                      type="email" 
                      value="john.smith@example.com"
                      className="w-full bg-gray-700 border border-gray-600 text-white p-3 rounded-lg"
                    />
                  </div>
                  
                  <div>
                    <label className="block text-white font-medium mb-2">Institution</label>
                    <input 
                      type="text" 
                      value="Research Institute of Molecular Science"
                      className="w-full bg-gray-700 border border-gray-600 text-white p-3 rounded-lg"
                    />
                  </div>
                  
                  <div className="pt-2">
                    <Button className="bg-blue-600 hover:bg-blue-700 text-white">
                      Save Changes
                    </Button>
                  </div>
                </div>
              </div>
              
              <div id="api" className="bg-gray-800/50 border border-gray-700 p-6 rounded-xl mb-6">
                <h2 className="text-white text-2xl font-bold mb-4">API Configuration</h2>
                <p className="text-gray-300 mb-4">Configure API settings for molecule generation and analysis.</p>
                
                <div className="space-y-4">
                  <div>
                    <label className="block text-white font-medium mb-2">API Key</label>
                    <div className="flex">
                      <input 
                        type="password" 
                        value="●●●●●●●●●●●●●●●●"
                        className="flex-1 bg-gray-700 border border-gray-600 text-white p-3 rounded-l-lg"
                        readOnly
                      />
                      <Button className="bg-gray-700 hover:bg-gray-600 text-white rounded-l-none border border-gray-600">
                        Show
                      </Button>
                    </div>
                  </div>
                  
                  <div>
                    <label className="block text-white font-medium mb-2">API Endpoint</label>
                    <input 
                      type="url" 
                      value="https://api.molecule-ai.example.com/v1"
                      className="w-full bg-gray-700 border border-gray-600 text-white p-3 rounded-lg"
                    />
                  </div>
                  
                  <div className="pt-2">
                    <Button className="bg-blue-600 hover:bg-blue-700 text-white">
                      Update API Settings
                    </Button>
                  </div>
                </div>
              </div>
              
              <div id="notifications" className="bg-gray-800/50 border border-gray-700 p-6 rounded-xl">
                <h2 className="text-white text-2xl font-bold mb-4">Notifications</h2>
                <p className="text-gray-300">Configure how you want to receive notifications.</p>
                
                <div className="mt-4 space-y-3">
                  <div className="flex items-center justify-between">
                    <label className="text-white">Email Notifications</label>
                    <div className="relative inline-block w-12 h-6 rounded-full bg-gray-700">
                      <input type="checkbox" className="sr-only" id="toggle-email" defaultChecked />
                      <span className="block absolute left-1 top-1 bg-blue-500 w-4 h-4 rounded-full transition-transform translate-x-6"></span>
                    </div>
                  </div>
                  
                  <div className="flex items-center justify-between">
                    <label className="text-white">Task Completion Alerts</label>
                    <div className="relative inline-block w-12 h-6 rounded-full bg-gray-700">
                      <input type="checkbox" className="sr-only" id="toggle-task" defaultChecked />
                      <span className="block absolute left-1 top-1 bg-blue-500 w-4 h-4 rounded-full transition-transform translate-x-6"></span>
                    </div>
                  </div>
                  
                  <div className="flex items-center justify-between">
                    <label className="text-white">System Updates</label>
                    <div className="relative inline-block w-12 h-6 rounded-full bg-gray-700">
                      <input type="checkbox" className="sr-only" id="toggle-updates" />
                      <span className="block absolute left-1 top-1 bg-gray-600 w-4 h-4 rounded-full transition-transform"></span>
                    </div>
                  </div>
                </div>
              </div>
            </div>
          </div>
        </main>
      </div>
      
      <Footer />
    </div>
  );
};

export default Settings;
